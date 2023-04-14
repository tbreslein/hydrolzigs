const std = @import("std");
const assert = std.debug.assert;

const hydrol_path = "src/hydrol.zig";

pub fn build(b: *std.build.Builder) void {
    const target = b.standardTargetOptions(.{});
    const mode = b.standardReleaseOptions();

    const sim_path = b.option(
        []const u8,
        "sim",
        "path to the .zig file for your simulation",
    ) orelse "";
    const opts = b.addOptions();
    opts.addOption([]const u8, "sim", sim_path);

    if (sim_path.len > 0) {
        var sim_path_splits = std.mem.splitBackwards(u8, sim_path, "/");
        const sim_file = sim_path_splits.next().?;
        assert(std.mem.eql(u8, sim_file[sim_file.len - 4 .. sim_file.len], ".zig"));
        const name = sim_file[0 .. sim_file.len - 4];

        const exe = b.addExecutable(name, sim_path);
        exe.setTarget(target);
        exe.setBuildMode(mode);
        exe.addPackagePath("hydrol", hydrol_path);
        // exe.install();

        const run_cmd = exe.run();
        run_cmd.step.dependOn(b.getInstallStep());
        if (b.args) |args| {
            run_cmd.addArgs(args);
        }
        const run_step = b.step(name, "Run this simulation");
        run_step.dependOn(&run_cmd.step);
    }

    const lib_tests = b.addTest(hydrol_path);
    lib_tests.setTarget(target);
    lib_tests.setBuildMode(mode);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&lib_tests.step);
}
