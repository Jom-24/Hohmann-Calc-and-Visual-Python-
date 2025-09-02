## Hohmann Transfer Calculator and Visualizer

A Python script by John Matthew Atok to calculate and visualize a Hohmann transfer orbit between two co-planar, circular orbits.

## Description
This project provides a tool to calculate the key parameters of a Hohmann transfer, including the required changes in velocity (delta-V) for each burn and the total time of flight. The script uses data from NASA fact sheets for various celestial bodies in our solar system, allowing for transfers around planets, moons, or the Sun. 

A key feature is the calculation of the Sphere of Influence (SOI) for the selected central body. The script will provide a warning if the target orbit exceeds this limit, as the model used for the calculations is only accurate within this sphere. The final output is a 2D plot of the initial orbit, the transfer ellipse, and the final orbit, generated using `matplotlib`.

## Features
- Calculates orbital velocities, energies, and the required delta-V for both burns.
- Computes the total time of flight for the transfer in seconds, minutes, hours, and days.
- Supports a wide range of central bodies, including the Sun, Earth, Moon, Mars, and more.
- Automatically calculates the Sphere of Influence (SOI) to ensure the accuracy of the model.
- Generates a clean and informative 2D plot of the entire maneuver for clear visualization.

## How to Use

Under the "inputs" comment, there are 3 variables that can be changed:
- central_body (list of bodies available are in the "celestial_bodies" library above the inputs comment)
- alt1 (this is the initial altitude of the orbit above the surface of the central body in km)
- alt2 (this is the final or desired orbit altitude above the surface of the central bodt in km)

There is no need to change the different values for each body (mu, radius, mass, etc) since the code automatically pulls the relevant information based on the central body that was set under "central_body".
