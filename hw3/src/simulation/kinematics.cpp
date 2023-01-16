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

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    // TODO
    // You need to return the solution (x) of the linear least squares system:
    //     i.e. find x which min(| Jacobian * x - target |)
    Eigen::VectorXd deltatheta = Jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(target);
    // Eigen::VectorXd deltatheta = (Jacobian.transpose() * Jacobian).inverse() * Jacobian.transpose()*target;
    
    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW3 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
                             acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;
    // TODO
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    size_t bone_num = 0;
    // HINT:
    // calculate number of bones need to move to perform IK, store in `bone_num`
    // a.k.a. how may bones from end_bone to its parent than to start_bone (include both side)
    acclaim::Bone* temp = end_bone;
    while (temp != start_bone && temp != nullptr) {
        //std::cout << temp->idx << std::endl;
        temp=temp->parent;
        bone_num++;
        
    }
    if (temp!=nullptr)
        bone_num++;
    // std::cout << "bone num: " << bone_num << std::endl;


    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
        if (desiredVector.norm() < epsilon) {
            break;
        }
        // HINT:
        // Calculate Jacobian, store in `Jacobian`
        temp = end_bone;
        for (long long i = 0; i < bone_num; i++) {

            Eigen::Vector4d r = end_bone->end_position - temp->start_position;

            if (temp->dofrx) 
                Jacobian.col(3 * i) = temp->rotation.matrix().col(0).cross3(r);
            if (temp->dofry) 
                Jacobian.col(3 * i + 1) = temp->rotation.matrix().col(1).cross3(r);

            if (temp->dofrz)
                Jacobian.col(3 * i + 2) = temp->rotation.matrix().col(2).cross3(r);

            temp = temp->parent;
        }

        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);
        // HINT:
        // Change `posture.bone_rotation` based on deltatheta
        temp = end_bone;
        for (long long i = 0; i < bone_num; i++) {



            posture.bone_rotations[temp->idx] =
                posture.bone_rotations[temp->idx] +
                util::toDegree(Eigen::Vector4d(deltatheta(3 * i) * temp->dofrx, deltatheta(3 * i + 1) * temp->dofry,
                                               deltatheta(3 * i + 2) * temp->dofrz, 0));


            temp = temp->parent;
        }

    }
    // TODO (Bonus)
    // Return IK is stable?
    // i.e. not swinging its hand in air
    if ((target_pos - end_bone->end_position).norm() < epsilon) {
        
        return true;
    } else {
        
        posture = original_posture;
        forwardSolver(posture, root_bone);
        return false;
    }

   
}
}  // namespace kinematics
