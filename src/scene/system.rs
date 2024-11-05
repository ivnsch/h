use std::f32::consts::PI;

use bevy::{color::palettes::css::GREEN, prelude::*};
use factorial::Factorial;
use num_complex::{Complex, Complex32, ComplexFloat};

pub fn add_scene(app: &mut App) {
    app.add_systems(Startup, add_dots);
}

pub fn add_dots(
    mut cmd: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let min = -5.;
    let max = 5.;
    let step = 0.05;

    let n = 1;
    let l = 0;
    let m = 0;

    let mut x = min;
    let mut y;
    let mut z;

    let mesh_handle = meshes.add(Cuboid { ..default() }.mesh());
    let material_handle = materials.add(StandardMaterial {
        base_color: GREEN.into(),
        ..default()
    });

    let mesh = Mesh3d(mesh_handle);
    let material = MeshMaterial3d(material_handle);

    let cube_scale = 0.05;
    let mut sphere_count = 0;
    while x <= max as f32 {
        x += step;
        y = min;
        while y <= max as f32 {
            y += step;
            z = min;
            while z <= max as f32 {
                z += step;

                let sc = to_spheric(Vec3::new(x, y, z));
                let val = psi_mod(sc, n, l, m).unwrap();

                if val > 0.1 {
                    let _ = cmd
                        .spawn((
                            mesh.clone(),
                            material.clone(),
                            Transform::from_translation(Vec3::new(x, y, z))
                                .with_scale(Vec3::new(cube_scale, cube_scale, cube_scale)),
                        ))
                        .id();
                }

                sphere_count += 1;
            }
        }
    }
    println!("finish render! objs: {}", sphere_count);
}

struct SphericCoords {
    rad: f32,
    theta: f32,
    phi: f32,
}

fn to_spheric(coords: Vec3) -> SphericCoords {
    let x = coords.x;
    let y = coords.y;
    let z = coords.z;
    let rad = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
    let theta = (y / x).atan();
    let phi = (z / rad).acos();

    SphericCoords { rad, theta, phi }
}

fn psi_mod(sc: SphericCoords, n: u8, l: u8, m: i8) -> Result<f32, String> {
    let psi = psi(sc, n, l, m)?;
    assert!(!psi.re.is_nan());
    assert!(!psi.norm_sqr().is_nan());
    Ok(psi.norm_sqr())
}

// https://en.wikipedia.org/wiki/Hydrogen_atom#Wavefunction
// under "normalized position wavefunctions"
fn psi(sc: SphericCoords, n: u8, l: u8, m: i8) -> Result<Complex32, String> {
    if n <= l {
        return Err(format!("n: {} <= l: {}", n, l));
    }

    // https://en.wikipedia.org/wiki/Bohr_radius#Reduced_Bohr_radius
    // let rbr = 5.29177210544e-11;
    // we can set the reduced bohr radius to 1, allowing to pass radius in the same
    let rbr = 1.;
    // floats casts here, for less verbosity
    let nf = n as f32;
    // let lf = l as f32;

    // left part under square root
    // these terms don't have any meaning other than internal grouping here
    let term1 = (2. / nf * rbr).powi(3);
    let term2 = (n - l - 1).factorial() as f32 / ((2 * n) * (n + l).factorial()) as f32;
    let term3 = (term1 * term2).sqrt();

    let p = (2. * sc.rad) / (nf * rbr);
    let term4 = std::f32::consts::E.powf(-p / 2.);
    let term5 = term4 * p.powi(l.into());
    let term6 = term3 * term5;

    // can do n - l - 1 because we checked n > l
    let term7 = laguerre_pol(n - l - 1, p)?;
    let term8 = spheric_harmonic(l, m, sc.theta, sc.phi)?;

    // println!("p: {}", p);
    // println!("term1: {}", term1);
    // println!("term2: {}", term2);
    // println!("term3: {}", term3);
    // println!("term4: {}", term4);
    // println!("term5: {}", term5);
    // println!("term6: {}", term6);
    // println!("term7: {}", term7);
    // println!("term8: {}", term8);
    assert!(!term1.is_nan());
    assert!(!term3.is_nan());
    assert!(!term4.is_nan());
    assert!(!term5.is_nan());
    assert!(!term6.is_nan());
    assert!(!term7.is_nan());
    assert!(!term8.re.is_nan());

    Ok(term6 * term7 * term8)
}

// https://en.wikipedia.org/wiki/Laguerre_polynomials#The_first_few_polynomials
fn laguerre_pol(n: u8, x: f32) -> Result<f32, String> {
    Ok(match n {
        0 => 1.,
        1 => -x + 1.,
        2 => 1. / 2. * (x.powi(2) - 4. * x + 2.),
        3 => 1. / 6. * (x.powi(3) + 9. * x.powi(2) - 18. * x + 6.),
        // TODO
        _ => return Err(format!("Not supported: n: {}", n)),
    })
}

// https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Complex_spherical_harmonics
fn spheric_harmonic(l: u8, m: i8, theta: f32, phi: f32) -> Result<Complex32, String> {
    // quick access
    let oh = 1. / 2.; // one half
    let e = std::f32::consts::E;
    let ex = phi * Complex::i(); // exponent (part)

    Ok(match (l, m) {
        (0, 0) => Complex::new(oh * (1. / PI).sqrt(), 0.),
        (1, -1) => oh * (3. / (2. * PI)).sqrt() * e.powc(-ex) * theta.sin(),
        (1, 0) => Complex::new(oh * (3. / PI).sqrt() * theta.cos(), 1.),
        (1, 1) => -oh * (3. / (2. * PI)).sqrt() * e.powc(ex) * theta.sin(),
        // TODO
        _ => return Err(format!("Not supported: l: {}, m: {}", l, m)),
    })
}

// // https://brilliant.org/wiki/spherical-harmonics/
// fn spheric_harmonic(l: u8, m: i8, theta: f32, phi: f32) -> Complex32 {
//     let mf = m as f32;
//     let lf = l as f32;

//     let term1 = (2. * lf + 1.) / (4. * PI);
//     let term2 = (l - m as u8).factorial() as f32 / (l + m as u8).factorial() as f32;
//     let term3 = (term1 * term2).sqrt();

//     // this has a derivative, so hardcoding the harmonics
//     let term4 = legrende_pol(l, m, theta.cos());

//     let exp = mf * phi * Complex::i();
//     let term5 = std::f32::consts::E.powc(exp);

//     term3 * term4 * term5
// }
