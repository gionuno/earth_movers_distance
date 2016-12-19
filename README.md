# earth_movers_distance

Calculates a dist matrix for the color histograms of a collection of images.

The images were taken from the Birds and Butterflies datasets: http://www-cvr.ai.uiuc.edu/ponce_grp/data/

First the histograms are made, using calc_hists.py.
This includes an interior point method in C++, using Armadillo, to solve the linear programming problem required to obtain said distance.

To compile the method, run 'setup.py build' in a terminal.

WARNING! Takes a VERY LONG TIME. But it's fun and seems like a good waste of time.
