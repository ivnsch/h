use crate::{
    camera_controller::{CameraController, CameraControllerPlugin},
    rotator::{Rotator, RotatorPlugin},
};
use bevy::prelude::*;

#[allow(dead_code)]
pub fn add_3d_space(app: &mut App) {
    app.add_plugins((DefaultPlugins, CameraControllerPlugin, RotatorPlugin))
        .add_systems(Startup, (setup_camera, setup_light));
}

fn setup_light(mut commands: Commands) {
    commands.insert_resource(AmbientLight {
        color: Color::WHITE,
        brightness: 1.0,
    });

    commands.spawn(DirectionalLight {
        illuminance: light_consts::lux::OVERCAST_DAY,
        shadows_enabled: true,
        ..default()
    });
}

fn setup_camera(mut commands: Commands) {
    commands.spawn((
        Camera3d { ..default() },
        Transform::from_xyz(0., 0., 8.0).looking_at(Vec3::ZERO, Vec3::Y),
        CameraController::default(),
        Rotator::default(),
    ));
}
