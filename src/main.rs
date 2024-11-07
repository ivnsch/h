mod camera_controller;
mod rotator;
mod scene;
mod system_3d;

use bevy::{app::App, math::Vec3};
use scene::system::add_scene;
use system_3d::add_3d_space;

// NOTE: this project tests a WGSL shader: syntax has to be compatible as far as possible

fn main() {
    let app = &mut App::new();
    add_3d_space(app);
    add_scene(app);
    app.run();

    // test_shader_code();
}

#[allow(unused)]
fn test_shader_code() {
    let x = 0.4;
    let y = 0.1;
    let z = 0.1;

    let n: u32 = 2;
    let l: u32 = 0;
    let m = 0;

    let point = Vec3::new(x, y, z);
    let spheric_coords = to_spheric_coords(point);
    let psi = psi(spheric_coords, n, l, m);
    let prob = prob(spheric_coords, n, l, m);

    println!("psi: {}, prob: {}", psi, prob);
}

fn sqrt(n: f32) -> f32 {
    n.sqrt()
}

fn pow(n: f32, e: f32) -> f32 {
    n.powf(e)
}

fn atan(n: f32) -> f32 {
    n.atan()
}

fn acos(n: f32) -> f32 {
    n.acos()
}

fn to_spheric_coords(coords: Vec3) -> SphericCoords {
    let rad = sqrt(pow(coords.x, 2.0) + pow(coords.y, 2.0) + pow(coords.z, 2.0));
    let theta = atan(coords.y / coords.x);
    let phi = acos(coords.z / rad);
    return SphericCoords { rad, theta, phi };
}

#[derive(Copy, Clone, Debug)]
pub struct SphericCoords {
    pub rad: f32,
    pub theta: f32,
    pub phi: f32,
}

pub fn prob(coords: SphericCoords, n: u32, l: u32, m: i32) -> f32 {
    let p = psi(coords, n, l, m);
    // we're returning just the real part so this is ok for now
    return pow(p, 2.0);
}

pub fn psi(coords: SphericCoords, n: u32, l: u32, m: i32) -> f32 {
    let rad = coords.rad;
    let theta = coords.theta;
    let phi = coords.phi;

    let e = 2.71828;
    // https://en.wikipedia.org/wiki/Bohr_radius#Reduced_Bohr_radius
    // let rbr = 5.29177210544e-11;
    // we can set the reduced bohr radius to 1, allowing to pass radius in the same
    let rbr = 1.;
    let nf = n as f32;

    // left part under square root
    // these terms don't have any meaning other than internal grouping here
    let term1 = pow(2. / nf * rbr, 3.0);
    let term2 = factorial(n - l - 1) as f32 / ((2 * n) * factorial(n + l)) as f32;
    let term3 = sqrt(term1 * term2);

    let p = (2. * rad) / (nf * rbr);
    let term4 = pow(e, -p / 2.0);
    let term5 = term4 * pow(p, l as f32);
    let term6 = term3 * term5;

    // can do n - l - 1 because we checked n > l
    let term7 = laguerre_pol(n - l - 1, p);
    let term8 = spheric_harmonic(l, m, theta, phi);

    // for now we'll just ignore the complex part
    return term6 * term7 * term8.real;
}

// https://en.wikipedia.org/wiki/Laguerre_polynomials#The_first_few_polynomials
fn laguerre_pol(n: u32, x: f32) -> f32 {
    match n {
        0 => {
            return 1.0;
        }
        1 => {
            return -x + 1.0;
        }
        2 => {
            return 1. / 2. * (pow(x, 2.0) - 4. * x + 2.);
        }
        3 => {
            return 1. / 6. * (pow(x, 3.0) + 9. * pow(x, 2.0) - 18. * x + 6.);
        }
        4 => {
            return 1. / 24. * (pow(x, 4.0) - 16. * pow(x, 3.0) + 72. * pow(x, 2.0) - 96. * x + 24.)
        }
        _ => {
            return 0.0;
        }
    }
}

// https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Complex_spherical_harmonics
// NOTE: don't improve syntax unless it's compatible with WGSL
fn spheric_harmonic(l: u32, m: i32, theta: f32, phi: f32) -> Complex {
    // quick access
    let oh = 1. / 2.; // one half

    let pi = 3.14159;

    let ex = Complex::new(0., -phi); // exponent (part)

    match l {
        0 => match m {
            0 => {
                return Complex::new(oh * sqrt(1. / pi), 0.);
            }
            _ => {
                return Complex::new(0., 0.);
            }
        },
        1 => {
            match m {
                -1 => {
                    return mul(e_to_i_pow(neg(ex)), oh * sqrt(3. / (2. * pi)) * sin(theta));
                }
                0 => {
                    let real = oh * sqrt(3. / pi) * cos(theta);
                    return Complex::new(real, 1.);
                }
                1 => {
                    return mul(e_to_i_pow(ex), -oh * sqrt(3. / (2. * pi)) * sin(theta));
                }
                _ => {
                    // TODO throw error?
                    return Complex::new(0., 0.);
                }
            }
        }
        2 => match m {
            -2 => {
                let ex = Complex::new(0., -2. * phi);
                return mul(
                    e_to_i_pow(ex),
                    1. / 4. * sqrt(15. / (2. * pi)) * pow(sin(theta), 2.),
                );
            }

            -1 => {
                let ex = Complex::new(0., -phi);
                return mul(
                    e_to_i_pow(ex),
                    oh * sqrt(15. / (2. * pi)) * sin(theta) * cos(theta),
                );
            }
            0 => {
                return Complex::new(
                    1. / 4. * sqrt(5. / pi) * (3. * pow(cos(theta), 2.) - 1.),
                    0.,
                );
            }
            1 => {
                let ex = Complex::new(0., phi);
                return mul(
                    e_to_i_pow(ex),
                    -oh * sqrt(15. / (2. * pi)) * sin(theta) * cos(theta),
                );
            }
            2 => {
                let ex = Complex::new(0., 2. * phi);
                return mul(
                    e_to_i_pow(ex),
                    1. / 4. * sqrt(15. / (2. * pi)) * pow(sin(theta), 2.),
                );
            }
            _ => {
                return Complex::new(0., 0.);
            }
        },

        _ => {
            return Complex::new(0., 0.);
        }
    }
}

struct Complex {
    real: f32,
    complex: f32,
}

impl Complex {
    fn new(real: f32, complex: f32) -> Complex {
        Complex { real, complex }
    }
}

fn e_to_i_pow(c: Complex) -> Complex {
    // e^ix = cos x + isin x,
    return Complex::new(cos(c.real), sin(c.complex));
}

fn neg(c: Complex) -> Complex {
    return Complex::new(-c.real, -c.complex);
}

fn mul(c: Complex, r: f32) -> Complex {
    return Complex::new(c.real * r, -c.complex * r);
}

fn sin(f: f32) -> f32 {
    f.sin()
}

fn cos(f: f32) -> f32 {
    f.cos()
}

fn factorial(n: u32) -> u32 {
    match n {
        0 => {
            return 1;
        }
        1 => {
            return 1;
        }
        2 => {
            return 2;
        }
        3 => {
            return 6;
        }
        4 => {
            return 24;
        }
        5 => {
            return 120;
        }
        // TODO error
        _ => {
            return 0;
        }
    }
}
