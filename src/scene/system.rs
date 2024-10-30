use bevy::{color::palettes::css::GREEN, prelude::*};

pub fn add_scene(app: &mut App) {
    app.add_systems(Startup, add_sphere);
}

pub fn add_sphere(
    mut cmd: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mesh_handle = meshes.add(Sphere { ..default() }.mesh().uv(32, 18));
    let material_handle = materials.add(StandardMaterial {
        base_color: GREEN.into(),
        ..default()
    });

    let sphere = Mesh3d(mesh_handle);
    let material = MeshMaterial3d(material_handle);

    let _ = cmd
        .spawn((
            sphere,
            material,
            Transform::from_xyz(0., 0., 0.).with_scale(Vec3::new(1., 1., 1.)),
        ))
        .id();
}
