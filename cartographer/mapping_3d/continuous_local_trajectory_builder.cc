#include <cartographer/mapping_3d/continuous_local_trajectory_builder.h>

namespace cartographer {
namespace mapping_3d {

ContinuousLocalTrajectoryBuilder::ContinuousLocalTrajectoryBuilder()
{

}
ContinuousLocalTrajectoryBuilder::ContinuousLocalTrajectoryBuilder(
    const proto::LocalTrajectoryBuilderOptions& options)
{

}
void ContinuousLocalTrajectoryBuilder::AddImuData(common::Time time, const Eigen::Vector3d& linear_acceleration,
                                             const Eigen::Vector3d& angular_velocity)
{

}
void ContinuousLocalTrajectoryBuilder::AddOdometerData(common::Time time,
                                                  const transform::Rigid3d& odometer_pose)
{

}
void ContinuousLocalTrajectoryBuilder::AddPlaneData(const common::Time time,
                                               const Eigen::Vector4d& coefficients)
{

}
}
}
