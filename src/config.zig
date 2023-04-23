// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports structs needed for configuring a hydrolzigs simulation.

/// Struct that configures a hydrolzigs simulation.
/// An instance of this struct should be available at comptime in order to help the compiler optimise as
/// aggressively as possible.
pub const Config = struct {
    /// Configures Mesh objects
    mesh: MeshConfig,

    /// Configures Physics objects
    physics: PhysicsConfig,

    /// Configures NumFlux objects
    numflux: NumFluxConfig,

    pub const MeshConfig = struct {
        /// Enumerates the different kinds of meshes available.
        pub const Type = enum {
            /// Regular cartesian mesh
            cartesian,
        };

        /// The kind of Mesh to construct (i.e. cartesian, logcylindrical, etc.)
        type: Type,

        /// The number of cells in the mesh, including ghost cells
        n: u32,

        /// xi coordinate at the inner or west edge of the computational area
        xi_in: f64,

        /// xi coordinate at the inner or east edge of the computational area
        xi_out: f64,

        /// Ensures that the different options are coherent. These checks are done at comptime and throw compile errors.
        fn validate(comptime self: MeshConfig) void {
            if (!(self.n > 5)) {
                @compileError("MeshConfig.n_all > 5 must hold!");
            }
            if (!(self.xi_in < self.xi_out)) {
                @compileError("MeshConfig.xi_in < MeshConfig.xi_out must hold!");
            }
        }
    };

    pub const PhysicsConfig = struct {
        /// Enumerates the different kinds of physics systems available.
        pub const Type = enum {
            /// adiabatic one-dimensional Euler equations
            euler1d_adiabatic,
            /// isothermal one-dimensional Euler equations
            euler1d_isothermal,
        };

        /// The type of physics to use
        type: Type,

        /// only used in adiabatic physics
        adiabatic_index: f64 = 1.0,
    };

    pub const NumFluxConfig = struct {
        /// Enumerates the different kinds of limiter functions available.
        pub const Limiter = enum {
            /// Min-Mod 1 limiter
            minmod,

            /// Superbee limiter
            superbee,

            /// Monocent limiter
            monocent,

            /// Van leer limiter
            vanleer,
        };

        /// The type of limiter function to use
        limiter_mode: Limiter = .vanleer,

        /// Generic parameter for those limiter functions that one (currently only Monocent)
        limiter_param: f64 = 1.2,
    };

    /// Ensures that the different options are coherent. These checks are done at comptime and throw compile errors.
    pub fn validate(comptime self: Config) void {
        self.mesh.validate();
    }
};
