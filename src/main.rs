mod camera_controller;
mod rotator;
mod scene;
mod system_3d;

use bevy::app::App;
use scene::system::add_scene;
use system_3d::add_3d_space;

fn main() {
    let app = &mut App::new();

    add_3d_space(app);

    add_scene(app);

    app.run();
}
