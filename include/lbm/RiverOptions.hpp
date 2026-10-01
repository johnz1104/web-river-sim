#pragma once

// River settings as text: "key=value" pairs that override RiverParams. The native
// preview (tools/river_preview.cpp) and the website (wasm/bindings.cpp) share this
// parser, so a setting tuned in one behaves the same in the other.
//
// Keys (numbers unless noted):
//   rows           lattice rows over the page height (integer, 8..1024)
//   cell           page units per lattice cell, 2..64 (sets rows; wins over rows)
//   outletCells    two-bank rivers: rows of outlet channel above the page, 0..400
//   re             Reynolds number Q / nu
//   speed          mean inflow speed, page units per second
//   lattice        mean inflow in lattice units (sets the Mach number and time step)
//   seed           random seed (unsigned integer)
//   model          collision model: bgk | trt | mrt
//   smag           Smagorinsky constant; 0 switches the closure off
//   sb, se         MRT bulk and even-ghost rates, 0 < s < 2
//   sq             MRT odd-ghost rate, 0 <= s < 2 (0 derives it from magic)
//   magic          TRT/MRT magic parameter
//   ramp, warmup   inflow ramp-up and hidden warm-up, seconds
//   inflowNoise    evolving inflow profile, 0 | 1
//   inflowAmplitude  how uneven it gets: std of mode m is amplitude / m (0.15)
//   inflowModes    how many cosine lanes across the inlet, 1..8 (3)
//   inflowFloor    the slowest a lane may get, as a fraction of the mean (0.3)
//   meander        transverse inflow meander, 0 | 1
//   stir           random stirring force, 0 | 1
//   stirIntensity  stirring strength relative to the inflow speed
//   kT             thermal noise temperature, lattice units (research only)
//   spacing        mean distance between dots, page units
//   diffusivity    dots' random-walk diffusivity, page units^2 per second
//   push           side push strength: velocity gained per second at its centre, in
//                  inlet speeds (0 = off; see PushParams in lbm/River.hpp)
//   pushX, pushY   the push's centre, page units
//   pushRadius     its Gaussian radius, page units > 0
//   pushAimX, pushAimY  the point it pushes toward, page units
//   pushPulse      pulsing: strength x (1 + pulse sin(2 pi t / period)), 0..1
//   pushPeriod     the pulse period, seconds > 0
//   fill           radius of a soft disc filled with dots when the dots start, so a
//                  swirl is populated at once (0 = off; see FillParams)
//   fillX, fillY   the disc's centre, page units
//   fillHold       seconds after the dots start during which the disc is topped back
//                  up to its first count (0 = fill once)
//   side           the designed river: the right edge takes water in below sideTop
//                  and lets it out above (0 | 1; see SideInflowParams)
//   sideTop        page y where the side inflow starts (outflow above it)
//   sideBottom     page y where it ends (a wall below; default: the page bottom)
//   sideSpeed      mean leftward speed of the side inflow, page units/s
//   sideAngle      degrees the side inflow tilts upward from leftward
//   sideMargin     page units the domain extends past the page's right edge
//   sideOutlet     the right edge above sideTop: 1 = outlet, 0 = wall
//   sideOutletBottom  with sideOutlet=1: page y where the right-edge outlet ends, with a
//                  wall from there down to sideTop (default: no wall)
//   topOutlet      designed river: 1 = outlet along the top (default), 0 = the top is a
//                  wall and water leaves only through the right edge (needs sideOutlet=1)
//   outletRight    page x where the top outlet ends (an invisible lid above the page
//                  right of it); default: no lid
//   boxX, boxY     top-left corner of a solid box in the water, page units (the page's
//                  control box sits on it; see ObstacleParams in River.hpp)
//   boxW, boxH     its size, page units (0 = no box)
//   cornerJet      width of a faster bottom-inflow stretch next to the left bank, page
//                  units (0 = off; see CornerJetParams in River.hpp)
//   cornerSpeed    its speed, page units/s
//   cornerAngle    degrees to the right of straight up (away from the bank), (-80, 80)
//   cornerEdge     page units over which it blends into the rest of the inflow
//   converge       aim the inflow at a focus point, 0 | 1 (see ConvergeParams in River.hpp)
//   convergeX, convergeY  the focus, page units (y < 0: above the page)
//   convergeMax    largest tilt toward it, degrees in (0, 89)
//   sideConverge   aim the side inflow at the focus too (instead of sideAngle), 0 | 1
//   depth          hidden bed depth as bottom drag, 0 | 1 (see DepthParams in River.hpp):
//                  layers for the inflow (1), the upper curve (2) and the outflow (3),
//                  and a shallow sill where layers 1 and 3 meet
//   bedDrag        drag rate at depth 1, per second (drag = bedDrag / depth^2)
//   depth1, depth2, depth3  the layers' relative depths, > 0
//   sillDepth      relative depth on the sill, > 0
//   sillX, sillY   the sill's tip; it runs from the right edge along y = sillY, page units
//   sillWidth      page units across the sill (0 = no sill)
//   curveBottom    page y of layer 2's lower edge (left of sillX)
//   stepWidth      page units over which the depth changes between layers

#include "lbm/River.hpp"

#include <string>

namespace lbm {

// Applies one setting. Returns "" on success, otherwise an error message; on error
// the parameters are unchanged.
std::string applyRiverOption(RiverParams& params, const std::string& key,
                             const std::string& value);

// Applies settings separated by spaces, commas or newlines, e.g. "speed=40 stir=1".
// Returns "" on success. On the first error it returns that error and leaves the
// parameters unchanged.
std::string applyRiverOptions(RiverParams& params, const std::string& list);

}  // namespace lbm
