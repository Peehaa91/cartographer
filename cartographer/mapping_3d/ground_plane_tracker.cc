#include "cartographer/mapping_3d/ground_plane_tracker.h"
#include "Eigen/Geometry"
#include "glog/logging.h"

namespace cartographer{
namespace mapping_3d{
GroundPlaneTracker::GroundPlaneTracker() :
		time_(common::Time::min()),
		last_coefficients_time_(common::Time::min()),
		orientation_(Eigen::Quaterniond(0,0,0,0))
{
}
GroundPlaneTracker::GroundPlaneTracker(Eigen::Quaterniond& quaternion, common::Time time ) :
		time_(time),
		last_coefficients_time_(common::Time::min()),
		orientation_(quaternion)
{

}

void GroundPlaneTracker::AddGroundPlaneObservation(const Eigen::Vector4d& coefficients)
{
	Eigen::Vector3d plane_normal(coefficients(0), coefficients(1), coefficients(2));
	LOG(INFO)<<"coeff: "<<plane_normal;
	Eigen::Vector3d z_axis(0, 0, 1);
	plane_normal.normalize();
	Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(plane_normal, z_axis);
	Eigen::Matrix3d rot_mat = Eigen::Matrix3d(quat);
	Eigen::Vector3d angles = rot_mat.eulerAngles(0, 1, 2);
	angles = Eigen::Vector3d(angles(0), angles(1), angles(2));
	Eigen::AngleAxisd rollAngle(angles(0), Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd pitchAngle(angles(1), Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd yawAngle(0, Eigen::Vector3d::UnitZ());
	Eigen::Quaterniond q = rollAngle * pitchAngle * yawAngle;
	Eigen::Matrix3d m = Eigen::Matrix3d(q);
	angles = m.eulerAngles(0, 1, 2);
	orientation_ = q;
	LOG(INFO)<<"roll: "<<angles(0)*180/M_PI<<" pitch: "<<angles(1)*180/M_PI<< " yaw: "<<angles(2)*180/M_PI;

	//LOG(INFO)<<"from calc:"<<orientation_.x()<<":"<<orientation_.y()<<":"<<orientation_.z()<<":"<<orientation_.w();

}

void GroundPlaneTracker::Advance(common::Time time)
{

}
}
}
