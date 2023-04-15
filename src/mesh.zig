const m_config = @import("config.zig");
const math = @import("std").math;

const n_gc_default: u32 = 2;

pub fn Mesh(comptime N: u32) type {
    return struct {
        comptime type: m_config.MeshConfig.Type = .cartesian,
        comptime n: u32 = N,
        comptime n_gc: u32 = n_gc_default,
        comptime n_comp: u32 = N - 2 * n_gc_default,
        comptime imin: u32 = 0,
        comptime ixi_in: u32 = n_gc_default,
        comptime ixi_out: u32 = N - n_gc_default - 1,
        xi_in: f64,
        xi_out: f64,
        dxi: f64,
        deta: f64,
        dphi: f64,
        xi_cent: @Vector(N, f64),
        xi_west: @Vector(N, f64),
        xi_east: @Vector(N, f64),
        xi_cent_inv: @Vector(N, f64),
        xi_west_inv: @Vector(N, f64),
        xi_east_inv: @Vector(N, f64),
        h_xi_cent: @Vector(N, f64),
        h_xi_west: @Vector(N, f64),
        h_xi_east: @Vector(N, f64),
        h_eta_cent: @Vector(N, f64),
        h_eta_west: @Vector(N, f64),
        h_eta_east: @Vector(N, f64),
        h_phi_cent: @Vector(N, f64),
        h_phi_west: @Vector(N, f64),
        h_phi_east: @Vector(N, f64),
        sqrt_g: @Vector(N, f64),
        line_xi: @Vector(N, f64),
        line_xi_inv: @Vector(N, f64),
        d_area_xi_deta_dphi_west: @Vector(N, f64),
        d_area_xi_deta_dphi_east: @Vector(N, f64),
        area_cell: @Vector(N, f64),
        area_west: @Vector(N, f64),
        area_east: @Vector(N, f64),
        volume: @Vector(N, f64),
        deta_dphi_d_volume: @Vector(N, f64),
        cell_width: @Vector(N, f64),
        cell_width_inv: @Vector(N, f64),
        cexe: @Vector(N, f64),
        cpxp: @Vector(N, f64),
        cxex: @Vector(N, f64),
        cpep: @Vector(N, f64),
        cxpx: @Vector(N, f64),
        cepe: @Vector(N, f64),
        minus_cexe: @Vector(N, f64),
        minus_cpxp: @Vector(N, f64),
        minus_cxex: @Vector(N, f64),
        minus_cpep: @Vector(N, f64),
        minus_cxpx: @Vector(N, f64),
        minus_cepe: @Vector(N, f64),
    };
}

pub fn init_mesh(comptime conf: m_config.Config) Mesh(conf.mesh.n) {
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
        .xi_in = conf.mesh.xi_in,
        .xi_out = conf.mesh.xi_out,
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
