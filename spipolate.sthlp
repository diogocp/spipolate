{smcl}
{vieweralsosee "[D] ipolate" "help ipolate"}{...}
{vieweralsosee "Sp introduction" "help sp_intro"}{...}
{viewerjumpto "Syntax" "spipolate##syntax"}{...}
{viewerjumpto "Description" "spipolate##description"}{...}

{p2colset 1 14 14 2}{...}
{p2col:{bf:spipolate} {hline 2}}Interpolation of spatial data{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 20 2}
{cmd:spipolate}
{varlist}
{cmd:using}
{it:{help filename}}
[{cmd:,}
{opt idw}[{cmd:(}{it:power}{cmd:)}]
{opt near:est}
{opt radius(#)}]


{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Methods}
{synopt :{opt idw}[{cmd:(}{it:power}{cmd:)}]}use the inverse distance weighting
	method with exponent {it:power}. This is the default method, and the
	default {it:power} is 2.{p_end}
{...}
{synopt :{opt near:est}}use the nearest neighbour method.{p_end}
{...}

{syntab: Options}
{synopt :{opt radius(#)}}specify the search radius; data points at a distance
	greater than this from the point being interpolated will be ignored.
{p_end}
{...}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
Only one method ({opt iwd} or {opt nearest}) may be specified.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{opt spipolate} creates the variables in {varlist} by spatial interpolation,
using data from {it:{help filename}}{cmd:.dta} (called the using dataset).

{pstd}
Both the dataset currently in memory and the using dataset are required to be
{bf:{help spset:[SP] spset}}. The variables in {varlist} must exist only in the
using dataset.
