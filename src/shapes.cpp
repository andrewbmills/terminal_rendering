#include "shapes.hpp"

Point rotate_point(Point p, const double R[9])
{
    Point p_out;
    p_out.x = R[0]*p.x + R[1]*p.y + R[2]*p.z;
    p_out.y = R[3]*p.x + R[4]*p.y + R[5]*p.z;
    p_out.z = R[6]*p.x + R[7]*p.y + R[8]*p.z;
    return p_out;
}

double dot_product(Point a, Point b) {
    double dot_product = (a.x*b.x + a.y*b.y + a.z*b.z);
    return dot_product;
}

Point normalize_vector(Point a) {
    Point b;
    double norm = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    b.x = a.x/norm;
    b.y = a.y/norm;
    b.z = a.z/norm;
    return b;
}

double normalized_dot_product(Point c, Point d) {
    Point a = normalize_vector(c);
    Point b = normalize_vector(d);
    double dot_product = (a.x*b.x + a.y*b.y + a.z*b.z);
    return dot_product;
}

Shape::Shape()
{
    set_orientation(0.0, 0.0, 0.0);
    m_screen.clear();
    BrightnessZBuffer b;
    m_screen.resize(m_pixel_height*m_pixel_width, b);
    set_orientation(0.0, 0.0, 0.0);
    set_light_source(0.0, -1.0, -1.0);
    return;
}

char Shape::convert_gray_to_ascii(const unsigned gray_scale_value)
{
    char ascii_gray = m_gray_to_ascii_12[gray_scale_value*(m_gray_to_ascii.size()-1)/255];
    return ascii_gray;
}

char Shape::convert_brightness_to_ascii(const double brightness)
{
    if (brightness == 0.0) return ' ';
    char ascii_gray = m_gray_to_ascii_12[round(brightness*(m_gray_to_ascii_12.size()-1))];
    return ascii_gray;
}

void Shape::set_orientation(const double roll, const double pitch, const double yaw)
{
    m_orientation.roll = roll;
    m_orientation.pitch = pitch;
    m_orientation.yaw = yaw;
    m_orientation.R[0] = cos(yaw)*cos(pitch);
    m_orientation.R[1] = cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll);
    m_orientation.R[2] = cos(yaw)*sin(pitch)*cos(roll) + sin(yaw)*sin(roll);
    m_orientation.R[3] = sin(yaw)*cos(pitch);
    m_orientation.R[4] = sin(yaw)*sin(pitch)*sin(roll) + cos(yaw)*cos(roll);
    m_orientation.R[5] = sin(yaw)*sin(pitch)*cos(roll) - cos(yaw)*sin(roll);
    m_orientation.R[6] = -sin(pitch);
    m_orientation.R[7] = cos(pitch)*sin(roll);
    m_orientation.R[8] = cos(pitch)*cos(roll);
    return;
}

void Shape::set_resolution(const unsigned height, const unsigned width)
{
    m_pixel_height = height;
    m_pixel_width = width;
    m_screen.clear();
    BrightnessZBuffer b;
    m_screen.resize(m_pixel_height*m_pixel_width, b);
    return;
}

void Shape::set_light_source(const double x, const double y, const double z)
{
    m_light_source.x = x;
    m_light_source.y = y;
    m_light_source.z = z;
    m_light_source = normalize_vector(m_light_source);
    return;
}

Pixel Shape::project_point_to_pixel(const Point p)
{
    Pixel pixel;
    double x_screen = p.x*m_K1/(m_observer_z - p.z);
    double y_screen = p.y*m_K1/(m_observer_z - p.z);
    if (std::abs(x_screen) < m_screen_width/2.0 && std::abs(y_screen) < m_screen_height) {
        pixel.u = m_pixel_width*(x_screen + m_screen_width/2.0)/m_screen_width;
        pixel.v = m_pixel_height*(m_screen_height/2.0 - y_screen)/m_screen_height;
    } else {
        pixel.u = m_pixel_width;
        pixel.v = m_pixel_height;
    }
    return pixel;
}

void Shape::render()
{
    BrightnessZBuffer b;
    m_screen.assign(m_screen.size(), b);
    for (int i=0; i<m_surface.size(); ++i) {
        Point p = rotate_point(m_surface[i], m_orientation.R);
        Point n = rotate_point(m_normals[i], m_orientation.R);
        Pixel pixel = project_point_to_pixel(p);
        if (pixel.u < m_pixel_width && pixel.v < m_pixel_height) {
            unsigned id = pixel.v*m_pixel_width + pixel.u;
            if (p.z > m_screen[id].z) {
                double brightness = std::abs(dot_product(n, m_light_source));
                m_screen[id].brightness = brightness;
                m_screen[id].z = p.z;
            }
        }
    }
    return;
}

void Shape::render_in_terminal()
{
    render();
    printf("\x1b[H");
    for (int i=0; i<m_screen.size(); ++i) {
        if (i % m_pixel_width == 0) std::cout << "\n";
        if (m_screen[i].brightness >= 0.0) std::cout << convert_brightness_to_ascii(m_screen[i].brightness);
        else std::cout << " ";
    }
    std::cout << std::endl;
    return;
}

Donut::Donut(const double t_R, const double t_r)
{
    m_R = t_R;
    m_r = t_r;
    m_K1 = m_screen_width*m_observer_z*3/(8*(m_r + m_R));
    compute_points_and_normals();
    return;
}

void Donut::compute_points_and_normals()
{
    for (double psi = 0.0; psi < 2*pi; psi = psi + m_psi_spacing) {
        for (double theta = 0.0; theta < 2*pi; theta = theta + m_theta_spacing) {
            Point p;
            p.x = cos(psi)*(m_R + m_r*cos(theta));
            p.y = m_r*sin(theta);
            p.z = sin(psi)*(m_R + m_r*cos(theta));
            m_surface.push_back(p);
            Point n;
            n.x = cos(psi)*cos(theta);
            n.y = sin(theta);
            n.z = sin(psi)*cos(theta);
            m_normals.push_back(n);
        }
    }
    return;
}

RectangularPrism::RectangularPrism(const double t_l, const double t_w, const double t_d)
{
    m_l = t_l;
    m_w = t_w;
    m_d = t_d;
    m_K1 = m_screen_width*m_observer_z*3/(8*1.414*std::max(m_l, std::max(m_w, m_d)));
    compute_points_and_normals();
    return;
}

void RectangularPrism::compute_points_and_normals()
{
    std::vector<Point> faces(6);
    faces[0].z = 0.5*m_d; // x-y plane
    faces[1].z = -0.5*m_d;
    faces[2].y = 0.5*m_l; // x-z plane
    faces[3].y = -0.5*m_l;
    faces[4].x = 0.5*m_w; // y-z plane
    faces[5].x = -0.5*m_w;

    std::vector<Point> face_normals(6);
    face_normals[0].z = 1.0;
    face_normals[1].z = -1.0;
    face_normals[2].y = 1.0;
    face_normals[3].y = -1.0;
    face_normals[4].x = 1.0;
    face_normals[5].x = -1.0;

    std::vector<double> dw;
    std::vector<double> dl;
    std::vector<double> dd;
    for (double width = -0.5*m_w; width <= 0.5*m_w; width += m_spacing) dw.push_back(width);
    for (double length = -0.5*m_l; length <= 0.5*m_l; length += m_spacing) dl.push_back(length);
    for (double depth = -0.5*m_d; depth <= 0.5*m_d; depth += m_spacing) dd.push_back(depth);
    

    for (int i=0; i<dw.size(); ++i) {
        for (int j=0; j<dl.size(); ++j) {
            faces[0].x = dw[i];
            faces[1].x = dw[i];
            faces[0].y = dl[j];
            faces[1].y = dl[j];
            m_surface.push_back(faces[0]);
            m_normals.push_back(face_normals[0]);
            m_surface.push_back(faces[1]);
            m_normals.push_back(face_normals[1]);
        }
    }
    for (int i=0; i<dw.size(); ++i) {
        for (int k=0; k<dd.size(); ++k) {
            faces[2].x = dw[i];
            faces[3].x = dw[i];
            faces[2].z = dd[k];
            faces[3].z = dd[k];
            m_surface.push_back(faces[2]);
            m_normals.push_back(face_normals[2]);
            m_surface.push_back(faces[3]);
            m_normals.push_back(face_normals[3]);
        }
    }
    for (int j=0; j<dl.size(); ++j) {        
        for (int k=0; k<dd.size(); ++k) {
            faces[4].y = dl[j];
            faces[5].y = dl[j];
            faces[4].z = dd[k];
            faces[5].z = dd[k];
            m_surface.push_back(faces[4]);
            m_normals.push_back(face_normals[4]);
            m_surface.push_back(faces[5]);
            m_normals.push_back(face_normals[5]);
        }
    }

    return;
}