use bevy::{input::mouse::MouseMotion, prelude::*, window::CursorGrabMode};

/// Based on Valorant's default sensitivity, not entirely sure why it is exactly 1.0 / 180.0,
/// but I'm guessing it is a misunderstanding between degrees/radians and then sticking with
/// it because it felt nice.
#[allow(unused)]
pub const RADIANS_PER_DOT: f32 = 1.0 / 180.0;

pub struct RotatorPlugin;

impl Plugin for RotatorPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(
            Update,
            run_camera_rotator, // , mouse_handler
        );
    }
}

#[allow(unused)]
#[derive(Component, Debug)]
pub struct Rotator {
    pub key_y: KeyCode,
    pub key_z: KeyCode,
    pub key_x: KeyCode,
    pub key_i: KeyCode,
    pub key_o: KeyCode,
    pub key_p: KeyCode,
    pub key_shift_left: KeyCode,
    pub key_shift_right: KeyCode,
}

impl Default for Rotator {
    fn default() -> Self {
        Self {
            key_y: KeyCode::KeyY,
            key_z: KeyCode::KeyZ,
            key_x: KeyCode::KeyX,
            key_i: KeyCode::KeyI,
            key_o: KeyCode::KeyO,
            key_p: KeyCode::KeyP,
            key_shift_left: KeyCode::ShiftLeft,
            key_shift_right: KeyCode::ShiftRight,
        }
    }
}

#[allow(clippy::too_many_arguments, unused)]
fn run_camera_rotator(
    key_input: Res<ButtonInput<KeyCode>>,
    mut camera: Query<(&mut Transform, &mut Rotator), With<Camera>>,
) {
    if let Ok((mut transform, mut rotator)) = camera.get_single_mut() {
        rotate(
            key_input,
            &rotator,
            &mut transform,
            rotator.key_x,
            rotator.key_y,
            rotator.key_z,
        );
    }
}

pub fn rotate(
    key_input: Res<ButtonInput<KeyCode>>,
    rotator: &Rotator,
    transform: &mut Transform,
    x_key: KeyCode,
    y_key: KeyCode,
    z_key: KeyCode,
) {
    // println!("rotating: {:?}", key_input);

    let mut rotation = 0.03;
    if key_input.pressed(rotator.key_shift_left) || key_input.pressed(rotator.key_shift_right) {
        rotation = -rotation;
    }

    if key_input.pressed(y_key) {
        transform.rotate_around(
            Vec3::ZERO,
            Quat::from_euler(EulerRot::XYZ, 0.0, rotation, 0.0),
        );
    }
    if key_input.pressed(z_key) {
        transform.rotate_around(
            Vec3::ZERO,
            Quat::from_euler(EulerRot::XYZ, 0.0, 0.0, rotation),
        );
    }
    if key_input.pressed(x_key) {
        transform.rotate_around(
            Vec3::ZERO,
            Quat::from_euler(EulerRot::XYZ, rotation, 0.0, 0.0),
        );
    }
}

/// updates cursor visibility and window focus for a given grab input
pub fn update_cursor_and_window_for_grab_input(
    windows: &mut Query<&mut Window>,
    mouse_events: &mut EventReader<MouseMotion>,
    input: &CursorGrabInput,
) {
    match input {
        CursorGrabInput::JustPressed => {
            for mut window in windows {
                if !window.focused {
                    continue;
                }
                window.cursor_options.grab_mode = CursorGrabMode::Locked;
                window.cursor_options.visible = false;
            }
        }
        CursorGrabInput::JustReleased => {
            for mut window in windows {
                window.cursor_options.grab_mode = CursorGrabMode::None;
                window.cursor_options.visible = true;
            }
            mouse_events.clear()
        }
    }
}

pub fn cursor_grab_update(
    mouse_button_input: Res<ButtonInput<MouseButton>>,
    button: MouseButton,
) -> Option<CursorGrabInput> {
    // Important: don't change ordering, sometimes pressed-release is delivered at the same time,
    // which we have to map to release, so we ask for release first
    if mouse_button_input.just_released(button) {
        return Some(CursorGrabInput::JustReleased);
    } else if mouse_button_input.just_pressed(button) {
        return Some(CursorGrabInput::JustPressed);
    }
    None
}

#[derive(Debug)]
pub enum CursorGrabInput {
    JustPressed,
    JustReleased,
}

#[derive(Debug)]
pub enum CursorGrabStatus {
    Active,
    Inactive,
}
