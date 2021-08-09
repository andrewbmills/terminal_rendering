#ifndef SHAPES_H
#define SHAPES_H

static const double pi = 3.1415926535;

#include <iostream>
#include <math.h>
#include <vector>

struct Orientation {
    double roll = 0.0;
    double pitch = 0.0;
    double yaw = 0.0;
    double R[9];
};

struct Pixel {
    unsigned u = 0;
    unsigned v = 0;
};

struct Point {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Line {
    Point origin;
    Point direction;
};

struct BrightnessZBuffer {
    double brightness = -1.0;
    double z = -1e6;
};

Line rotate_line(Line line, Orientation rotation);
double dot_product(Point a, Point b);
Point normalize_vector(Point a);
double normalized_dot_product(Point a, Point b);

class Shape
{
    public:
        const std::string m_gray_to_ascii = "$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/|()1{}[]?-_+~<>i!lI;:,\"^`'. ";
        const std::string m_gray_to_ascii_12 = ".,-~:;=!*#$@";
        Orientation m_orientation;
        unsigned m_pixel_width = 100; // window width in pixels for rendering
        unsigned m_pixel_height = 100; // window height in pixels for rendering
        double m_observer_z = 5.0; // distance from shape origin to observer eye
        double m_screen_z = 50.0; // distance from shape origin to screen
        double m_screen_width = 100.0;
        double m_screen_height = 100.0;
        double m_K1;
        Point m_light_source;
        std::vector<BrightnessZBuffer> m_screen;
        std::vector<Point> m_surface;
        std::vector<Point> m_normals;
        Shape();
        char convert_gray_to_ascii(const unsigned gray_scale_value);
        char convert_brightness_to_ascii(const double brightness);
        void set_orientation(const double roll, const double pitch, const double yaw);
        void set_resolution(const unsigned width, const unsigned height);
        void set_light_source(const double x, const double y, const double z);
        Pixel project_point_to_pixel(const Point p);
        virtual void render();
        virtual void render_in_terminal();
        // virtual void computer_points_and_normals() = 0;
};

class Donut: public Shape
{
    private:
        double m_R; // radius from donut center to loop center
        double m_r; // outer loop radius
        double m_psi_spacing = 0.1;
        double m_theta_spacing = 0.08;
        // double m_psi_spacing = 0.175;
        // double m_theta_spacing = 0.175;
        void compute_points_and_normals();
    public:
        Donut(const double t_R, const double t_r);
};

class RectangularPrism: public Shape
{
    private:
        double m_l; // length
        double m_w; // width
        double m_d; // depth
        double m_spacing = 0.04;
        void compute_points_and_normals();
    public:
        RectangularPrism(const double t_l, const double t_w, const double t_d);
};

#endif