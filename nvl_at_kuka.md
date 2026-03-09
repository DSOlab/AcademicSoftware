## General

The (lidar-based) NVL component is installed at `f400u04@10.112.112.44`. It resides under the folder `kuka_ws`.
`kuka_ws` will be hence named `ROOT_DIR`.

## Launch Scripts

| Path | Dir | Script | Launches | Published (topics/TF) | Consumes (Topics/TF) | Configs |
|------|-----|--------|----------|-----------------------|----------------------|---------|
|`soprano_slam_navigation`| `launch` | `launch_robot.launch.py` | `rsp.launch.py`, `sensor.launch.py`, `twist_mux`, `plc_bridge`| `twist_mux` → `/cmd_vel` (from mux) | `/cmd_vel_joy`, `/cmd_vel_tracker`, `/cmd_vel_nav` (inputs to mux); plus whatever `sensor.launch.py` | CLI args: `use_sim_time`, rosbag; RViz/URDF/etc. via `rsp.launch.py` |
|`nvl_component`          | `launch` | `sensor.launch.py` | `lidar1`, `lidar2`, `scan_merger`, `rf2o`| `/scan_merged` (merged LaserScan); `/odom_rf2o` (nav_msgs/Odometry); **TF**: `odom → base_link`; (in live mode also the two raw scans from urg drivers) | `/lidar_front_right/scan`, `/lidar_rear_left/scan` (from drivers or rosbag); `/clock` when `use_sim_time:=true`; **TF** lookups to `base_link` & lidar frames for merging | `params_ether_fr.yaml`, `params_ether_rl.yaml` (urg_node2); (merger params inline in launch); **CLI** args: `rosbag`, `use_sim_time` |
|`soprano_slam_navigation`| `launch` | `slam_navigation.launch.py` | `slam_toolbox`, `navigation_launch.py`| `slam_toolbox` → `/map`, `/map_metadata`, `/pose` (PoseWithCovarianceStamped), **TF** (map linkage if enabled); `Nav2` → planners/controllers topics (e.g., /cmd_vel_nav), costmaps, action servers | Laser scan (as configured, e.g., `/scan_merged`), odometry/TF from NVL (`odom → base_link`), `/clock` when `use_sim_time:=true` | `slam_toolbox.yaml` (SLAM), `nav2_params.yaml` (Nav2), **CLI** args: `use_sim_time`, `slam_params_file`, `params_file`|
|`soprano_slam_navigation`| `launch` | `rsp.launch.py` | `robot_state_publisher` | **TF** tree from URDF (`/tf`, `/tf_static` as applicable); sets robot_description params | | `description/robot.urdf.xacro` (includes: `robot_base.xacro`, `lidar.xacro`, `d456.urdf.xacro`, `ros2_control.xacro`); CLI arg: `use_sim_time`|


## Installation

Go to `ROOT_DIR` and type:

```
colcon build --symlink-install
source install/setup.bash
```

## Configuration Files

### Geometry (Mounting Points) & TFs

Static configuration (i.e. LiDAR mounting points and static **TF**s w.r.t. payload and `base_link`) are published via `rsp.launch.py` and consumed by all other scripts. Hence, it is important to make sure that **xacro/urdf** info are correctly published. These files are located at `/src/soprano_slam_navigation/description/`.

Both `slam_navigation.launch.py` and `sensor.launch.py` use these info.

### Lidar IPs and (static) Configuration

these have to be set in the config files:
`${ROOT_DIR}/src/nvl_component/config/params_ether_fr.yaml` for the front right Lidar and 
`${ROOT_DIR}/src/nvl_component/config/params_ether_rl.yaml` for the rear left respectively.

Current (default) values for the Lidars are:

* 192.168.10.61:10940 for the front-right, and
* 192.168.10.62:10940 for the rear-left.

Note that their frames must be named `lidar_front_right_frame` and `lidar_rear_left_frame`. 
Below is an example for the front-right Lidar:

```
f400u04@ubuntu:~/kuka_ws$ cat src/nvl_component/config/params_ether_fr.yaml 
urg_node2:
  ros__parameters:
    ip_address: '192.168.10.61'
    ip_port: 10940
    # Match URDF frame name
    frame_id: 'lidar_front_right_frame'
    
    calibrate_time: false
    synchronize_time: false
    publish_intensity: false
    publish_multiecho: false
    error_limit: 4
    error_reset_period: 5.0
    diagnostics_tolerance: 0.05
    diagnostics_window_time: 5.0
    time_offset: 0.0
    angle_min: -2.35619   # -135 deg
    angle_max: 2.35619    # 135 deg
    skip: 0
    cluster: 1
```

## Testing With Rosbags

## Launch Scripts

### `src/nvl_component/launch/sensor.launch.py`

Runs either live LiDARs or a rosbag (mutually exclusive). 
- If you pass a **rosbag path** (`ros2 launch ... rosbag:=/path/to/bag`), it **does not** start the two `urg_node2` LiDAR drivers and instead **plays the bag** (after a 2.5 s delay) with `/clock` at 50 Hz and only the two scan topics (`/lidar_front_right/scan`, `/lidar_rear_left/scan`).
- If you **don’t** pass a rosbag path, it starts **two lifecycle `urg_node2` nodes** (front/rear), configures & activates them, and remaps their outputs to the two scan topics.
- **Merges the two LaserScans** with `laser_merger2` into **`/scan_merged`** (using the inline parameters in the file).
- **Runs RF2O** on **`/scan_merged`** to publish **`/odom_rf2o`** and the **odom→base_link TF** (since `publish_tf: True`).
- Honors **sim time** (if enabled) for the merger and RF2O.

**Command-line arguments (name → default → what they do)**

- `rosbag` → `''`
  If **non-empty**, triggers **bag mode**: suppresses the live LiDAR nodes and runs `ros2 bag play <rosbag> --clock 50` for the two scan topics. If empty, runs live drivers.
- `use_sim_time` → `'false'`
  When `true`, sets `use_sim_time` on the merger and RF2O so they consume `/clock` (useful in bag mode). 
- `urg_config` → `<pkg>/config/params_ether_fr.yaml`
  YAML for the **front** `urg_node2`. Loaded and passed when in live-LiDAR mode. 
- `urg_config2` → `<pkg>/config/params_ether_rl.yaml`
  YAML for the **rear** `urg_node2`. Loaded and passed in live-LiDAR mode. 
- `auto_start` → `'true'`
  Declared but effectively unused (the lifecycle handlers don’t check it; they’re gated only by rosbag/not). Harmless. 
- `node_name` → `'urg_node2_front'`
  Name for the **front** LiDAR node. Used only in live-LiDAR mode. 
- `node_name2` → `'urg_node2_rear'`
  Name for the **rear** LiDAR node. Used only in live-LiDAR mode. 
- `scan_topic_name` → `'/lidar_front_right/scan'`
  Remap target for the **front** scan output (live-LiDAR mode). The merger also subscribes to this. 
- `scan_topic_name2` → `'/lidar_rear_left/scan'`
  Remap target for the **rear** scan output (live-LiDAR mode). The merger also subscribes to this. 

**Example invocations**

- **Live sensors** (default):
  `ros2 launch nvl_component sensor.launch.py`
- **Play a rosbag** (uses sim time & suppresses drivers):
  `ros2 launch nvl_component sensor.launch.py rosbag:=/data/my_run/ use_sim_time:=true`
  
# Real Examples

`ros2 launch soprano_slam_navigation launch_robot.launch.py use_sim_time:=true rosbag:=/home/xanthos/EvaluationData/test3_bag`
`ros2 launch soprano_slam_navigation slam_navigation.launch.py use_sim_time:=true`

## Upload to KUKA

- `src/nvl_component/launch/sensor.launch.py`
- `src/soprano_slam_navigation/launch/launch_robot.launch.py`
