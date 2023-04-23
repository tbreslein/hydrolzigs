// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

const std = @import("std");
const assert = std.debug.assert;

const hydrol_path = "src/hydrol.zig";

pub fn build(b: *std.build.Builder) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const hydrol_module = b.createModule(.{
        .source_file = .{ .path = hydrol_path },
    });

    const sim_path = b.option(
        []const u8,
        "sim",
        "path to the .zig file for your simulation",
    ) orelse "";
    const opts = b.addOptions();
    opts.addOption([]const u8, "sim", sim_path);

    if (sim_path.len > 0) {
        // Extract the file name without the ending, for example extract "template" from "examples/template.zig"
        // by splitting it into path chunks, go to the last chunk, and strip out the file ending.

        var sim_path_splits = std.mem.splitBackwards(u8, sim_path, std.fs.path.sep_str);
        const sim_file = sim_path_splits.next().?;
        assert(std.mem.eql(u8, sim_file[sim_file.len - 4 .. sim_file.len], ".zig"));
        const name = sim_file[0 .. sim_file.len - 4];

        const exe = b.addExecutable(.{
            .name = name,
            .root_source_file = .{ .path = sim_path },
            .target = target,
            .optimize = optimize,
        });
        exe.addModule("hydrol", hydrol_module);

        const run_cmd = b.addRunArtifact(exe);
        run_cmd.step.dependOn(b.getInstallStep());
        if (b.args) |args| {
            run_cmd.addArgs(args);
        }
        const run_step = b.step(name, "Run this simulation");
        run_step.dependOn(&run_cmd.step);

        b.default_step = run_step;
    }

    const lib_tests = b.addTest(.{
        .name = "hydrol_unit_test",
        .root_source_file = .{ .path = hydrol_path },
        .target = target,
        .optimize = optimize,
    });

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&lib_tests.step);
}
