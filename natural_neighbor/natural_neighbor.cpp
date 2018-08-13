/* Copyright (c) 2018 Diogo Pereira
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include "stplugin.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation;
typedef K::FT Coord_type;
typedef K::Point_2 Point;
typedef std::map<Point, Coord_type, K::Less_xy_2> Coord_map;


std::pair<std::vector<Coord_map>, ST_retcode> read_known_data(char *stata_matrix);


STDLL stata_call(int argc, char *argv[]) {
    if (argc != 1) return 197; // invalid syntax
    if (SF_nvars() < 3) return 102; // too few variables specified

    // Read known data from Stata
    std::pair<std::vector<Coord_map>, ST_retcode> data = read_known_data(argv[0]);
    if (data.second) return data.second;

    const std::vector<Coord_map> &value_functions = data.first;

    // Create Delaunay triangulations
    // We need one per interpolated variable because the variables may have different missing rows
    std::vector<Delaunay_triangulation> triangulations(value_functions.size());
    for (int v = 0; v < value_functions.size(); v++) {
        for (auto p : value_functions[v]) {
            triangulations[v].insert(p.first);
        }
    }

    // Interpolation main loop
    for (ST_int j = SF_in1(), sf_in2 = SF_in2(); j <= sf_in2; j++) {
        if (SF_ifobs(j)) {
            // Read coordinates of the interpolation point
            ST_double x, y;
            if (ST_retcode rc = SF_vdata(1, j, &x)) return rc;
            if (ST_retcode rc = SF_vdata(2, j, &y)) return rc;

            Point p = Point(x, y);

            // For each point, interpolate each variable v
            for (int v = 0; v < value_functions.size(); v++) {
                // Calculate the natural neighbor coordinates
                std::vector<std::pair<Point, Coord_type> > coords;
                auto result = CGAL::natural_neighbor_coordinates_2(triangulations[v], p, std::back_inserter(coords));

                if (result.third) {
                    // Interpolate
                    auto norm = result.second;
                    Coord_type res = CGAL::linear_interpolation(coords.begin(), coords.end(), norm,
                                                                CGAL::Data_access<Coord_map>(value_functions[v]));

                    // Send result to Stata
                    SF_vstore(v + 3, j, static_cast<ST_double>(res));
                }
            }
        }
    }

    return 0;
}


// Read data from a Stata matrix, where the first two columns of the matrix are the coordinates,
// and the remaining columns are the variables to be interpolated
std::pair<std::vector<Coord_map>, ST_retcode> read_known_data(char *stata_matrix) {
    ST_double x, y, val;
    ST_int num_vars = SF_col(stata_matrix) - 2; // The first two columns are the coordinates

    if (num_vars <= 0) return std::make_pair(std::vector<Coord_map>(), 503);

    std::vector<Coord_map> value_functions((unsigned) num_vars);

    for (ST_int i = 1; i <= SF_row(stata_matrix); i++) {
        if (ST_retcode rc = SF_mat_el(stata_matrix, i, 1, &x)) return std::make_pair(value_functions, rc);
        if (ST_retcode rc = SF_mat_el(stata_matrix, i, 2, &y)) return std::make_pair(value_functions, rc);

        if (SF_is_missing(x) || SF_is_missing(y)) continue;

        for (ST_int j = 0; j < num_vars; j++) {
            if (ST_retcode rc = SF_mat_el(stata_matrix, i, j + 3, &val)) return std::make_pair(value_functions, rc);

            if (!SF_is_missing(val)) {
                value_functions[j][Point(x, y)] = val;
            }
        }
    }

    return std::make_pair(value_functions, 0);
}
