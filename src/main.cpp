#include <thread>
#include <chrono>
#include <stdlib.h>
#include "shapes.hpp"

int main()
{
    // Donut s(2.0, 1.0);
    RectangularPrism s(2.0, 2.0, 2.0);
    s.set_resolution(40, 80);
    double roll = pi/4.0;
    double yaw = pi/4.0;
    double pitch = 0.0;
    s.set_orientation(roll, pitch, yaw);
    s.render_in_terminal();
    while (true) {
        s.render_in_terminal();
        roll = roll + 5.0*(pi/180.0);
        yaw = yaw - 2.5*(pi/180.0);
        pitch = pitch + 3.0*(pi/180.0);
        s.set_orientation(roll, pitch, yaw);
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }

    return 0;
}