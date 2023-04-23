// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the Mesh struct for building the grids a simulation runs on.

const Config = @import("config.zig").Config;
const math = @import("std").math;

const n_gc_default: u32 = 2;

/// A grid of cells, identified by coordinates xi, that a simulation runs on.
pub fn Mesh(comptime c: Config) type {
    const dxi = comptime switch (c.mesh.type) {
        .cartesian => (c.mesh.xi_out - c.mesh.xi_in) / @as(f64, c.mesh.n - 4),
    };
    const deta = comptime switch (c.mesh.type) {
        .cartesian => 1.0,
    };
    const dphi = comptime switch (c.mesh.type) {
        .cartesian => 1.0,
    };

    const xi_west = blk: {
        var vec = @splat(c.mesh.n, @as(f64, 0.0));
        comptime switch (c.mesh.type) {
            .cartesian => {
                inline for (0..c.mesh.n) |i| {
                    vec[i] = @mulAdd(f64, dxi, @as(i32, i) - 2, c.mesh.xi_in);
                }
            },
        };
        break :blk vec;
    };
    const xi_cent = xi_west + @splat(c.mesh.n, 0.5 * dxi);
    const xi_east = xi_west + @splat(c.mesh.n, dxi);

    const h_xi_cent = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_xi_west = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_xi_east = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_eta_cent = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_eta_west = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_eta_east = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_phi_cent = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_phi_west = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };
    const h_phi_east = comptime switch (c.mesh.type) {
        .cartesian => @splat(c.mesh.n, @as(f64, 1.0)),
    };

    const sqrt_g = h_xi_cent * h_eta_cent * h_phi_cent;

    const line_xi = h_xi_cent * @splat(c.mesh.n, dxi);
    const line_xi_inv = @splat(c.mesh.n, @as(f64, 1.0)) / line_xi;

    const d_area_xi_deta_dphi_west = h_eta_west * h_phi_west;
    const d_area_xi_deta_dphi_east = h_eta_east * h_phi_east;

    const area_cell = (xi_east - xi_west) * (xi_east * xi_west);
    const area_west = comptime d_area_xi_deta_dphi_west * @splat(c.mesh.n, @as(f64, deta * dphi));
    const area_east = d_area_xi_deta_dphi_east * @splat(c.mesh.n, @as(f64, deta * dphi));

    const volume = sqrt_g * @splat(c.mesh.n, dxi * deta * dphi);
    const deta_dphi_d_volume = @splat(c.mesh.n, @as(f64, deta * dphi)) / (@splat(c.mesh.n, @as(f64, math.f64_epsilon)) + volume);

    const cell_width = xi_east - xi_west;
    const cell_width_inv = @splat(c.mesh.n, @as(f64, 1.0)) / cell_width;

    const cexe = @splat(c.mesh.n, @as(f64, 0.5)) * (h_phi_east + h_phi_west) * (h_eta_east - h_eta_west) * deta_dphi_d_volume;
    const cpxp = @splat(c.mesh.n, @as(f64, 0.5)) * (h_eta_east + h_eta_west) * (h_phi_east - h_phi_west) * deta_dphi_d_volume;
    const cxex = @splat(c.mesh.n, @as(f64, 1.0));
    const cpep = @splat(c.mesh.n, @as(f64, 1.0));
    const cxpx = @splat(c.mesh.n, @as(f64, 1.0));
    const cepe = @splat(c.mesh.n, @as(f64, 1.0));

    return struct {
        /// The type of mesh; defaults to cartesian
        comptime type: Config.MeshConfig.Type = .cartesian,

        /// The number of ghost cells per edge; defaults to 2
        comptime n_gc: u32 = n_gc_default,

        /// The number of cells in the computational area; defaults to n - 2*n_gc
        comptime n_comp: u32 = c.mesh.n - 2 * n_gc_default,

        /// The first index in the computational area
        comptime ixi_in: u32 = n_gc_default,

        /// The last index in the computational area
        comptime ixi_out: u32 = c.mesh.n - n_gc_default - 1,

        /// The differential along the xi coordinatedoc
        comptime dxi: f64 = dxi,

        /// The differential along the eta coordinate
        comptime deta: f64 = deta,

        /// The differential along the Phi coordinate
        comptime dphi: f64 = dphi,

        /// Vector of xi coordinates in the centre of each cell
        comptime xi_cent: @Vector(c.mesh.n, f64) = xi_cent,

        /// Vector of xi coordinates at the west border of each cell
        comptime xi_west: @Vector(c.mesh.n, f64) = xi_west,

        /// Vector of xi coordinates at the east border of each cell
        comptime xi_east: @Vector(c.mesh.n, f64) = xi_east,

        /// Inverse of xi_cent
        comptime xi_cent_inv: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 1.0)) / xi_cent,

        /// Inverse of xi_west
        comptime xi_west_inv: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 1.0)) / xi_west,

        /// Inverse of xi_east
        comptime xi_east_inv: @Vector(c.mesh.n, f64) = @splat(c.mesh.n, @as(f64, 1.0)) / xi_east,

        /// Geometric scale along the xi coordinate at the centre of each cell
        comptime h_xi_cent: @Vector(c.mesh.n, f64) = h_xi_cent,

        /// Geometric scale along the xi coordinate at the west border of each cell
        comptime h_xi_west: @Vector(c.mesh.n, f64) = h_xi_west,

        /// Geometric scale along the xi coordinate at the east border of each cell
        comptime h_xi_east: @Vector(c.mesh.n, f64) = h_xi_east,

        /// Geometric scale along the eta coordinate at the centre of each cell
        comptime h_eta_cent: @Vector(c.mesh.n, f64) = h_eta_cent,

        /// Geometric scale along the eta coordinate at the west border of each cell
        comptime h_eta_west: @Vector(c.mesh.n, f64) = h_eta_west,

        /// Geometric scale along the eta coordinate at the east border of each cell
        comptime h_eta_east: @Vector(c.mesh.n, f64) = h_eta_east,

        /// Geometric scale along the Phi coordinate at the centre of each cell
        comptime h_phi_cent: @Vector(c.mesh.n, f64) = h_phi_cent,

        /// Geometric scale along the Phi coordinate at the west border of each cell
        comptime h_phi_west: @Vector(c.mesh.n, f64) = h_phi_west,

        /// Geometric scale along the Phi coordinate at the east border of each cell
        comptime h_phi_east: @Vector(c.mesh.n, f64) = h_phi_east,

        /// Square root of the metric factor g
        comptime sqrt_g: @Vector(c.mesh.n, f64) = sqrt_g,

        /// Line element along xi
        comptime line_xi: @Vector(c.mesh.n, f64) = line_xi,

        /// Inverse of line_xi
        comptime line_xi_inv: @Vector(c.mesh.n, f64) = line_xi_inv,

        /// Shorthand for the west cell face area perpendicular to xi divided by deta and dphi
        comptime d_area_xi_deta_dphi_west: @Vector(c.mesh.n, f64) = d_area_xi_deta_dphi_west,

        /// Shorthand for the east cell face area perpendicular to xi divided by deta and dphi
        comptime d_area_xi_deta_dphi_east: @Vector(c.mesh.n, f64) = d_area_xi_deta_dphi_east,

        /// Surface area of each cell in the eta-Phi plane
        comptime area_cell: @Vector(c.mesh.n, f64) = area_cell,

        /// Surface area for the west cell face area perpendicular to xi
        comptime area_west: @Vector(c.mesh.n, f64) = area_west,

        /// Surface area for the east cell face area perpendicular to xi
        comptime area_east: @Vector(c.mesh.n, f64) = area_east,

        /// Volume of each cell
        comptime volume: @Vector(c.mesh.n, f64) = volume,

        /// Shorthand for: deta * dphi / volume
        comptime deta_dphi_d_volume: @Vector(c.mesh.n, f64) = deta_dphi_d_volume,

        /// Distance between west and east ends of each cell along xi
        comptime cell_width: @Vector(c.mesh.n, f64) = cell_width,

        /// Inverse of cell_width
        comptime cell_width_inv: @Vector(c.mesh.n, f64) = cell_width_inv,

        /// Commutator coefficient for eta-xi-eta
        comptime cexe: @Vector(c.mesh.n, f64) = cexe,

        /// Commutator coefficient for Phi-xi-Phi
        comptime cpxp: @Vector(c.mesh.n, f64) = cpxp,

        /// Commutator coefficient for xi-eta-xi
        comptime cxex: @Vector(c.mesh.n, f64) = cxex,

        /// Commutator coefficient for Phi-eta-Phi
        comptime cpep: @Vector(c.mesh.n, f64) = cpep,

        /// Commutator coefficient for xi-Phi-xi
        comptime cxpx: @Vector(c.mesh.n, f64) = cxpx,

        /// Commutator coefficient for eta-Phi-eta
        comptime cepe: @Vector(c.mesh.n, f64) = cepe,

        /// Shorthand for -cexe
        comptime minus_cexe: @Vector(c.mesh.n, f64) = -cexe,

        /// Shorthand for -cpxp
        comptime minus_cpxp: @Vector(c.mesh.n, f64) = -cpxp,

        /// Shorthand for -cxex
        comptime minus_cxex: @Vector(c.mesh.n, f64) = -cxex,

        /// Shorthand for -cpep
        comptime minus_cpep: @Vector(c.mesh.n, f64) = -cpep,

        /// Shorthand for -cxpx
        comptime minus_cxpx: @Vector(c.mesh.n, f64) = -cxpx,

        /// Shorthand for -cepe
        comptime minus_cepe: @Vector(c.mesh.n, f64) = -cepe,
    };
}
