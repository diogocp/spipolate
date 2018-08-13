program define spipolate
	version 15
	syntax newvarlist using/, [radius(numlist min=1 max=1 >=0 missingokay) ///
		NEARest idw IDWpower(numlist min=1 max=1) NATural]

	if "`idwpower'" != "" local idw "idw"

	local nopts : word count `idw' `nearest'
	if `nopts' > 1 {
		display as error "must specify one interpolation method only"
		exit 198
	}

	mata: my_spipolate = spipolate()
	if "`nearest'" != "" mata: my_spipolate.method = "nearest"
	if "`natural'" != "" mata: my_spipolate.method = "natural"
	if "`radius'" != "" mata: my_spipolate.radius = `radius'
	if "`idwpower'" != "" mata: my_spipolate.idwpower = `idwpower'

	// Import the using data into Mata
	preserve
	use "`using'", nolabel
	confirm numeric variable `varlist', exact
	mata: my_spipolate.load_known_data("`varlist'")
	restore, preserve

	// Generate new variables to store the interpolation results
	foreach var in `varlist' {
		quietly generate `var' = .
	}

	mata: my_spipolate.interpolate()

	// All OK, cancel restore
	restore, not
end


version 15
program spipolate_natural_neighbor, plugin

set matastrict on
mata:

class spipolate {
	public:
		void new()
		void interpolate()
		void load_known_data()

		string scalar method
		real scalar radius, idwpower

	private:
		real colvector calc_distances()
		real rowvector idw()
		real rowvector nearest()
		void natural()

		string scalar coordsys
		string scalar varlist
		real matrix known_coords, known_data
}

void spipolate::new() {
	this.method = "idw"
	this.radius = .
	this.idwpower = 2
}

void spipolate::load_known_data(string scalar varlist) {
	string scalar cx, cy

	stata("quietly spset")
	cx = st_global("r(sp_cx)")
	cy = st_global("r(sp_cy)")
	this.coordsys = st_global("r(sp_coord_sys)")

	// Ignore points with missing coordinates
	stata("quietly drop if missing(" + cx + ") | missing(" + cy + ")")

	this.known_coords = st_data(., (cx, cy))
	this.known_data = st_data(., varlist)
	this.varlist = varlist
}

void spipolate::interpolate() {
	real matrix points, output
	real scalar i

	st_view(points, ., ("_CX", "_CY"))
	st_view(output, ., (this.varlist))

	if(method == "idw") {
		for (i = 1; i <= rows(points); i++) {
			output[i, .] = idw(points[i, .])
		}
	} else if(method == "nearest") {
		for (i = 1; i <= rows(points); i++) {
			output[i, .] = nearest(points[i, .])
		}
	} else if(method == "natural") {
		natural()
	} else {
		_error(3300, "method not implemented")
	}
}

real rowvector spipolate::idw(real rowvector x) {
	real colvector d

	d = calc_distances(x)

	// Check if there are any zero distances
	if (min(d) == 0) {
		// If yes, then we want every data point at zero distance to be
		// given the same weight, and ignore all other points. This can
		// be done with a logical not, which flips zeros to ones, and
		// changes all non-zero values (including missing) to zero.
		d = !d
	}

	// Ignore distances greater than the maximum radius
	if(radius != .) {
		d = d :* (d :<= radius)
	}

	return(colmeans(known_data, 1 :/ (d :^ idwpower)))
}

real rowvector spipolate::nearest(real rowvector x) {
	real colvector d
	real scalar min_d

	d = calc_distances(x)
	min_d = min(d)

	// If the nearest neighbor is outside the radius, return missing
	if(min_d > radius) {
		return(J(1, cols(known_data), .))
	}

	// There might be more than one point at the minimum distance,
	// so return an average of those
	return(colmeans(known_data, d :== min(d)))
}

void spipolate::natural() {
	string scalar st_known_data

	st_known_data = st_tempname()
	st_matrix(st_known_data, (known_coords, known_data))

	stata("plugin call spipolate_natural_neighbor latitude longitude " +
		varlist + ", " + st_known_data)
}

/* Returns a colvector with the distances from x to each of the
 * points in known_coords.
 */
real colvector spipolate::calc_distances(real rowvector x) {
	real colvector d

	if (coordsys == "latlong") {
		d = haversine_distance(x, known_coords)
	} else {
		_error(3300, "coordsys not implemented")
	}

	return(d)
}

/* Calculates the great-circle distance between points A and points B,
 * using the haversine formula.
 */
real colvector haversine_distance(real matrix A, real matrix B) {
	real matrix a, b, h

	// Convert from degrees to radians
	a = A * (pi() / 180)
	b = B * (pi() / 180)

	// If one of the arguments is a rowvector and the other is a matrix,
	// expand the rowvector to conform with the matrix.
	if(rows(a) == 1 & rows(b) > 1) {
		a = J(rows(b), 1, a)
	} else if(rows(a) > 1 & rows(b) == 1) {
		b = J(rows(a), 1, b)
	}

	// haversine(a-b)
	h = sin((a-b)/2) :^ 2

	// 12742: mean diameter of the Earth (km)
	return(12742 * asin(sqrt(h[,1] + cos(a[,1]) :* cos(b[,1]) :* h[,2])))
}

/* colmeans() is essentially the same as the standard mean(), except that
 * whereas mean(X, w) ignores rows of X with missing values, colmeans only
 * ignores the missing values themselves, not their entire rows.
 */
real rowvector colmeans(real matrix X, real colvector w) {
	real rowvector result
	real scalar i

	result = J(1, cols(X), .)

	for (i = 1; i <= cols(X); i++) {
		result[i] = mean(X[., i], w)
	}

	return(result)
}

end
