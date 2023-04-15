//! Exports the Mesh struct for building the grids a simulation runs on.

const m_config = @import("config.zig");
const math = @import("std").math;

const n_gc_default: u32 = 2;

/// A grid of cells, identified by coordinates xi, that a simulation runs on.
pub fn Mesh(comptime N: u32) type {
    return struct {
        /// The type of mesh; defaults to cartesian
        comptime type: m_config.MeshConfig.Type = .cartesian,

        /// The number of ghost cells per edge; defaults to 2
        comptime n_gc: u32 = n_gc_default,

        /// The number of cells in the computational area; defaults to n - 2*n_gc
        comptime n_comp: u32 = N - 2 * n_gc_default,

        /// The first index in the computational area
        comptime ixi_in: u32 = n_gc_default,

        /// The last index in the computational area
        comptime ixi_out: u32 = N - n_gc_default - 1,

        /// The differential along the xi coordinatedoc
        dxi: f64,

        /// The differential along the eta coordinate
        deta: f64,

        /// The differential along the Phi coordinate
        dphi: f64,

        /// Vector of xi coordinates in the centre of each cell
        xi_cent: @Vector(N, f64),

        /// Vector of xi coordinates at the west border of each cell
        xi_west: @Vector(N, f64),

        /// Vector of xi coordinates at the east border of each cell
        xi_east: @Vector(N, f64),

        /// Inverse of xi_cent
        xi_cent_inv: @Vector(N, f64),

        /// Inverse of xi_west
        xi_west_inv: @Vector(N, f64),

        /// Inverse of xi_east
        xi_east_inv: @Vector(N, f64),

        /// Geometric scale along the xi coordinate at the centre of each cell
        h_xi_cent: @Vector(N, f64),

        /// Geometric scale along the xi coordinate at the west border of each cell
        h_xi_west: @Vector(N, f64),

        /// Geometric scale along the xi coordinate at the east border of each cell
        h_xi_east: @Vector(N, f64),

        /// Geometric scale along the eta coordinate at the centre of each cell
        h_eta_cent: @Vector(N, f64),

        /// Geometric scale along the eta coordinate at the west border of each cell
        h_eta_west: @Vector(N, f64),

        /// Geometric scale along the eta coordinate at the east border of each cell
        h_eta_east: @Vector(N, f64),

        /// Geometric scale along the Phi coordinate at the centre of each cell
        h_phi_cent: @Vector(N, f64),

        /// Geometric scale along the Phi coordinate at the west border of each cell
        h_phi_west: @Vector(N, f64),

        /// Geometric scale along the Phi coordinate at the east border of each cell
        h_phi_east: @Vector(N, f64),

        /// Square root of the metric factor g
        sqrt_g: @Vector(N, f64),

        /// Line element along xi
        line_xi: @Vector(N, f64),

        /// Inverse of line_xi
        line_xi_inv: @Vector(N, f64),

        /// Shorthand for the west cell face area perpendicular to xi divided by deta and dphi
        d_area_xi_deta_dphi_west: @Vector(N, f64),

        /// Shorthand for the east cell face area perpendicular to xi divided by deta and dphi
        d_area_xi_deta_dphi_east: @Vector(N, f64),

        /// Surface area of each cell in the eta-Phi plane
        area_cell: @Vector(N, f64),

        /// Surface area for the west cell face area perpendicular to xi
        area_west: @Vector(N, f64),

        /// Surface area for the east cell face area perpendicular to xi
        area_east: @Vector(N, f64),

        /// Volume of each cell
        volume: @Vector(N, f64),

        /// Shorthand for: deta * dphi / volume
        deta_dphi_d_volume: @Vector(N, f64),

        /// Distance between west and east ends of each cell along xi
        cell_width: @Vector(N, f64),

        /// Inverse of cell_width
        cell_width_inv: @Vector(N, f64),

        /// Commutator coefficient for eta-xi-eta
        cexe: @Vector(N, f64),

        /// Commutator coefficient for Phi-xi-Phi
        cpxp: @Vector(N, f64),

        /// Commutator coefficient for xi-eta-xi
        cxex: @Vector(N, f64),

        /// Commutator coefficient for Phi-eta-Phi
        cpep: @Vector(N, f64),

        /// Commutator coefficient for xi-Phi-xi
        cxpx: @Vector(N, f64),

        /// Commutator coefficient for eta-Phi-eta
        cepe: @Vector(N, f64),

        /// Shorthand for -cexe
        minus_cexe: @Vector(N, f64),

        /// Shorthand for -cpxp
        minus_cpxp: @Vector(N, f64),

        /// Shorthand for -cxex
        minus_cxex: @Vector(N, f64),

        /// Shorthand for -cpep
        minus_cpep: @Vector(N, f64),

        /// Shorthand for -cxpx
        minus_cxpx: @Vector(N, f64),

        /// Shorthand for -ccepe
        minus_cepe: @Vector(N, f64),
    };
}

/// Initialises a Mesh instance
pub fn initMesh(comptime conf: m_config.Config) Mesh(conf.mesh.n) {
    const dxi = comptime switch (conf.mesh.type) {
        .cartesian => (conf.mesh.xi_out - conf.mesh.xi_in) / @as(f64, conf.mesh.n - 4),
    };
    const deta = comptime switch (conf.mesh.type) {
        .cartesian => 1.0,
    };
    const dphi = comptime switch (conf.mesh.type) {
        .cartesian => 1.0,
    };

    const xi_west = blk: {
        var vec = @splat(conf.mesh.n, @as(f64, 0.0));
        comptime switch (conf.mesh.type) {
            .cartesian => {
                // TODO: make this a ranged for loop with zig-git
                var i = 0;
                while (i < conf.mesh.n) : (i += 1) {
                    vec[i] = @mulAdd(f64, dxi, i - 2, conf.mesh.xi_in);
                }
            },
        };
        break :blk vec;
    };
    const xi_cent = xi_west + @splat(conf.mesh.n, 0.5 * dxi);
    const xi_east = xi_west + @splat(conf.mesh.n, dxi);

    const h_xi_cent = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_xi_west = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_xi_east = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_eta_cent = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_eta_west = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_eta_east = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_phi_cent = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_phi_west = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };
    const h_phi_east = comptime switch (conf.mesh.type) {
        .cartesian => @splat(conf.mesh.n, @as(f64, 0.0)),
    };

    const sqrt_g = h_xi_cent * h_eta_cent * h_phi_cent;

    const line_xi = h_xi_cent * @splat(conf.mesh.n, dxi);
    const line_xi_inv = @splat(conf.mesh.n, @as(f64, 1.0)) / line_xi;

    const d_area_xi_deta_dphi_west = h_eta_west * h_phi_west;
    const d_area_xi_deta_dphi_east = h_eta_east * h_phi_east;

    const area_cell = (xi_east - xi_west) * (xi_east * xi_west);
    const area_west = comptime d_area_xi_deta_dphi_west * @splat(conf.mesh.n, @as(f64, deta * dphi));
    const area_east = d_area_xi_deta_dphi_east * @splat(conf.mesh.n, @as(f64, deta * dphi));

    const volume = sqrt_g * @splat(conf.mesh.n, dxi * deta * dphi);
    const deta_dphi_d_volume = @divExact(@splat(conf.mesh.n, @as(f64, deta * dphi)), @splat(conf.mesh.n, @as(f64, math.f64_epsilon)) + volume);

    const cell_width = xi_east - xi_west;
    const cell_width_inv = @splat(conf.mesh.n, @as(f64, 1.0)) / cell_width;

    const cexe = @splat(conf.mesh.n, @as(f64, 0.5)) * (h_phi_east + h_phi_west) * (h_eta_east - h_eta_west) * deta_dphi_d_volume;
    const cpxp = @splat(conf.mesh.n, @as(f64, 0.5)) * (h_eta_east + h_eta_west) * (h_phi_east - h_phi_west) * deta_dphi_d_volume;
    const cxex = @splat(conf.mesh.n, @as(f64, 1.0));
    const cpep = @splat(conf.mesh.n, @as(f64, 1.0));
    const cxpx = @splat(conf.mesh.n, @as(f64, 1.0));
    const cepe = @splat(conf.mesh.n, @as(f64, 1.0));

    return Mesh(conf.mesh.n){
        .type = conf.mesh.type,
        .dxi = dxi,
        .deta = deta,
        .dphi = dphi,
        .xi_cent = xi_cent,
        .xi_west = xi_west,
        .xi_east = xi_east,
        .xi_cent_inv = @splat(conf.mesh.n, @as(f64, 1.0)) / xi_cent,
        .xi_west_inv = @splat(conf.mesh.n, @as(f64, 1.0)) / xi_west,
        .xi_east_inv = @splat(conf.mesh.n, @as(f64, 1.0)) / xi_east,
        .h_xi_cent = h_xi_cent,
        .h_xi_west = h_xi_west,
        .h_xi_east = h_xi_east,
        .h_eta_cent = h_eta_cent,
        .h_eta_west = h_eta_west,
        .h_eta_east = h_eta_east,
        .h_phi_cent = h_phi_cent,
        .h_phi_west = h_phi_west,
        .h_phi_east = h_phi_east,
        .sqrt_g = sqrt_g,
        .line_xi = line_xi,
        .line_xi_inv = line_xi_inv,
        .d_area_xi_deta_dphi_west = d_area_xi_deta_dphi_west,
        .d_area_xi_deta_dphi_east = d_area_xi_deta_dphi_east,
        .area_cell = area_cell,
        .area_west = area_west,
        .area_east = area_east,
        .volume = volume,
        .deta_dphi_d_volume = deta_dphi_d_volume,
        .cell_width = cell_width,
        .cell_width_inv = cell_width_inv,
        .cexe = cexe,
        .cpxp = cpxp,
        .cxex = cxex,
        .cpep = cpep,
        .cxpx = cpxp,
        .cepe = cepe,
        .minus_cexe = -cexe,
        .minus_cpxp = -cpxp,
        .minus_cxex = -cxex,
        .minus_cpep = -cpep,
        .minus_cxpx = -cxpx,
        .minus_cepe = -cepe,
    };
}
