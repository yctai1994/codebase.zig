wksp: []f64,
n: usize,
ncv: usize,
cnvgd: bool,
sum: f64,
eps: f64,
lastval: f64,
lasteps: f64,

const Error = error{OutOfWorkSpace};

const Self: type = @This();

fn init(allocator: mem.Allocator, nmax: usize, epss: f64) !*Self {
    const self: *Self = try allocator.create(Self);
    errdefer allocator.destroy(self);

    self.wksp = try allocator.alloc(f64, nmax);

    self.n = 0;
    self.ncv = 0;
    self.cnvgd = false;
    self.sum = 0.0;
    self.eps = epss;
    self.lastval = 0.0;

    return self;
}

fn deinit(self: *Self, allocator: mem.Allocator) void {
    allocator.free(self.wksp);
    allocator.destroy(self);
}

fn next(self: *Self, term: f64) Error!f64 {
    if (self.n + 1 > self.wksp.len) return error.OutOfWorkSpace;

    if (self.n == 0) { // Initialize
        self.wksp[self.n] = term;
        self.sum = 0.5 * term; // Return first estimate.
        self.n += 1;
    } else {
        var tmp: f64 = undefined;
        var dum: f64 = undefined;

        tmp = self.wksp[0];
        self.wksp[0] = term;

        // Update saved quantities by van Wijngaarden’s algorithm.
        for (1..self.n) |j| {
            dum = self.wksp[j];
            self.wksp[j] = 0.5 * (self.wksp[j - 1] + tmp);
            tmp = dum;
        }
        self.wksp[self.n] = 0.5 * (self.wksp[self.n - 1] + tmp);

        if (@abs(self.wksp[self.n]) <= @abs(self.wksp[self.n - 1])) {
            // Favorable to increase p, and the table becomes longer.
            self.sum += 0.5 * self.wksp[self.n];
            self.n += 1;
        } else {
            // Favorable to increase n, the table doesn’t become longer.
            self.sum += self.wksp[self.n];
        }
    }

    self.lasteps = @abs(self.sum - self.lastval);
    if (self.lasteps <= self.eps) self.ncv += 1;
    if (self.ncv >= 2) self.cnvgd = true;

    self.lastval = self.sum;
    return self.lastval;
}

test "test" {
    const page = testing.allocator;

    var obj: *Self = try Self.init(page, 50, 1e-16);
    defer obj.deinit(page);

    var tot: f64 = undefined;
    var ind: usize = 1;

    while (!obj.cnvgd) : (ind += 1) {
        tot = try obj.next(
            delta_potential_alter_k(ind, 5.0, 2.5, 1e-16),
        );
    }

    debug.print("{d}\n{any}\n", .{ tot, obj });
}

fn delta_potential_k(k: usize, d: f64, z: f64) f64 {
    const kd: f64 = @as(f64, @floatFromInt(k)) * d;
    return (2.0 / kd) - (1.0 / (kd + z)) - (1.0 / (kd + d - z));
}

fn delta_potential_alter_k(k: usize, d: f64, z: f64, tol: f64) f64 {
    var factor: usize = 1;
    var last_tot: f64 = 0.0;
    var last_val: f64 = math.inf(f64);

    while (tol <= @abs(last_val)) : (factor <<= 1) {
        last_val = @as(f64, @floatFromInt(factor)) * delta_potential_k(factor * k, d, z);
        last_tot += last_val;
    }

    return if ((k & 1) == 1) last_tot else -last_tot;
}

const std = @import("std");
const mem = std.mem;
const math = std.math;
const debug = std.debug;
const testing = std.testing;
