#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>
#include "acclaim/bone.h"
#include "util/helper.h"

namespace kinematics {

void computeBone2Global(const acclaim::Posture& posture, acclaim::Bone* bone) { 

    if (bone == nullptr) return;

    
    acclaim::Bone* parent = bone->parent;

    Eigen::Quaterniond local_rotation = util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);
    local_rotation.normalize();
    Eigen::Affine3d rot_local;
    rot_local = local_rotation.toRotationMatrix();

    bone->rotation = parent->rotation * bone->rot_parent_current * rot_local;
    bone->start_position = parent->end_position;
    bone->end_position = bone->start_position + bone->rotation * (bone->length * bone->dir);

    computeBone2Global(posture, bone->child);
    computeBone2Global(posture, bone->sibling);

}
void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO
    // This function will be called with bone == root bone of the skeleton
    // You should set these variables:
    //     bone->start_position = Eigen::Vector4d::Zero();
    //     bone->end_position = Eigen::Vector4d::Zero();
    //     bone->rotation = Eigen::Matrix4d::Zero();
    // The sample above just set everything to zero

    // traverse from root
    Eigen::Quaterniond root_rotation = util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);
    root_rotation.normalize();
    Eigen::Affine3d rot_root;
    rot_root = root_rotation.toRotationMatrix();

    bone->rotation = rot_root;
    bone->start_position = posture.bone_translations[bone->idx];
    bone->end_position = posture.bone_translations[bone->idx];
    
    computeBone2Global(posture, bone->child);
    computeBone2Global(posture, bone->sibling);


}

std::vector<acclaim::Posture> timeWarper(const std::vector<acclaim::Posture>& postures, int keyframe_old,
                                         int keyframe_new) {
    int total_frames = static_cast<int>(postures.size());
    int total_bones = static_cast<int>(postures[0].bone_rotations.size());
    std::vector<acclaim::Posture> new_postures = postures;
    for (int i = 0; i < total_frames; ++i) {
        for (int j = 0; j < total_bones; ++j) {
            // TODO
            // You should set these variables:
            //     new_postures[i].bone_translations[j] = postures[i].bone_translations[j];
            //     new_postures[i].bone_rotations[j] = postures[i].bone_rotations[j];
            // The sample above just change nothing

            if (i <= keyframe_new) {

                int frame1 = i * keyframe_old / keyframe_new;
                int frame2 = frame1 + 1;


                float ratio2 = (float)i * keyframe_old / keyframe_new - frame1;
                float ratio1 = 1 - ratio2;

                //std::cout << frame1 << " " << frame2 << " " << ratio1 << " " << ratio2 << std::endl;

                new_postures[i].bone_translations[j] =
                    postures[frame1].bone_translations[j] * ratio1 + postures[frame2].bone_translations[j] * ratio2;


                Eigen::Quaterniond rotation1 = util::rotateDegreeZYX(postures[frame1].bone_rotations[j]);
                rotation1.normalize();
                Eigen::Quaterniond rotation2 = util::rotateDegreeZYX(postures[frame2].bone_rotations[j]);
                rotation2.normalize();

                Eigen::Quaterniond interpolation = rotation1.slerp(ratio2, rotation2);
                interpolation.normalize();

                Eigen::Vector3d temp = interpolation.toRotationMatrix().eulerAngles(2, 1, 0);

                new_postures[i].bone_rotations[j] = util::toDegree(Eigen::Vector4d(temp[2], temp[1], temp[0], 0));

                
            
            } else {
                
                int frame1 = keyframe_old +
                             (i - keyframe_new) * (total_frames - keyframe_old - 1) / (total_frames - keyframe_new - 1);
                int frame2 = std::min(frame1 + 1, total_frames - 1);


                float ratio2 =
                    keyframe_old +
                    (float)(i - keyframe_new) * (total_frames - keyframe_old - 1) / (total_frames - keyframe_new - 1) -
                    frame1;
                float ratio1 = 1 - ratio2;

                //std::cout << frame1 << " " << frame2 << " " << ratio1 << " " << ratio2 << std::endl;

                new_postures[i].bone_translations[j] =
                    postures[frame1].bone_translations[j] * ratio1 + postures[frame2].bone_translations[j] * ratio2;


                Eigen::Quaterniond rotation1 = util::rotateDegreeZYX(postures[frame1].bone_rotations[j]);
                rotation1.normalize();
                Eigen::Quaterniond rotation2 = util::rotateDegreeZYX(postures[frame2].bone_rotations[j]);
                rotation2.normalize();

                Eigen::Quaterniond interpolation = rotation1.slerp(ratio2, rotation2);
                interpolation.normalize();

                Eigen::Vector3d temp = interpolation.toRotationMatrix().eulerAngles(2, 1, 0);

                new_postures[i].bone_rotations[j] = util::toDegree(Eigen::Vector4d(temp[2], temp[1], temp[0], 0));
            
            }

            
            
            

            // std::cout << interpolation.vec().size() << std::endl;
            
            /*new_postures[i].bone_rotations[j] = Eigen::Vector4d(interpolation.vec(), 0);*/
        }
    }
    return new_postures;
}
}  // namespace kinematics
