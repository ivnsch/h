use bevy::{color::palettes::css::GREEN, prelude::*};

use crate::{prob, to_spheric_coords};

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

    let n = 2;
    let l = 1;
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

                let vec = Vec3::new(x, y, z);
                let sc = to_spheric_coords(vec);
                let val = prob(sc, n, l, m);

                // some thresholds to "see something" with reasonable performance
                // for n=1, l=0, m=0
                // if val > 0.1 {
                // for 2s the probability for a given point is very low!
                // for n=2, l=0, m=0
                // if val > 0.002 {
                // for n=2, l=1, m=0
                if val > 0.026 {
                    let _ = cmd
                        .spawn((
                            mesh.clone(),
                            material.clone(),
                            Transform::from_translation(Vec3::new(x, y, z))
                                .with_scale(Vec3::new(cube_scale, cube_scale, cube_scale)),
                        ))
                        .id();
                    sphere_count += 1;
                }
            }
        }
    }
    println!("finish render! objs: {}", sphere_count);
}

#[cfg(test)]
mod test {
    use bevy::math::Vec3;

    use crate::{psi, to_spheric_coords};

    #[test]
    fn test() {
        let x = 0.1;
        let y = 0.;
        let z = 0.;

        let n: u32 = 2;
        let l: u32 = 0;
        let m = 0;

        let sc = to_spheric_coords(Vec3::new(x, y, z));
        println!("spheric: {:?}", sc);
        let val = psi(sc, n, l, m);

        assert_eq!(val, 0.08538426);
    }
}
