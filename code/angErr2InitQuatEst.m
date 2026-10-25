function q_est0 = angErr2InitQuatEst(eul_err, true_q_0)
    true_euler0 = quat2eul(true_q_0', 'ZYX');
    est_euler0 = true_euler0' +  eul_err;
    q_est0 = eul2quat(est_euler0', 'ZYX')';
end