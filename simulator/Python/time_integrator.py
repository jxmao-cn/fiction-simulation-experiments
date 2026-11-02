import copy
from cmath import inf

import numpy as np
import numpy.linalg as LA
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve

import Python_slovement.InertiaEnergy as InertiaEnergy
import Python_slovement.MassSpringEnergy as MassSpringEnergy
import Python_slovement.GravityEnergy as GravityEnergy
import Python_slovement.UniversalGravitationEnergy as UniversalGravitationEnergy
import Python_slovement.BarrierEnergy as BarrierEnergy
import Python_slovement.FrictionEnergy as FrictionEnergy
import Python_slovement.SpringEnergy as SpringEnergy

def step_forward(x, e, v, m, l2, k, n, o, contact_area, mu, is_DBC, DBC, DBC_v, DBC_limit, DBC_stiff, h, tol, *, stop_mode="step", grad_tol_rel=1e-7,
                 grad_scale=1.0, actual_step_tol=1e-12, actual_step_count=2,
                 max_iter=None, report=None):
    # 默认step保留原参考程序；其余模式用于显式诊断。
    if stop_mode not in ("step", "grad", "actual", "combined"):
        raise ValueError("invalid stop_mode")
    if not all(np.isfinite(a) and a > 0 for a in (h, tol, grad_tol_rel, grad_scale, actual_step_tol)) or not isinstance(actual_step_count, (int, np.integer)) or actual_step_count < 1 or (max_iter is not None and (not isinstance(max_iter, (int, np.integer)) or max_iter < 0)):
        raise ValueError("invalid stopping configuration")
    small_count = 0
    reason = ""
    force = float("nan")
    x_tilde = x + v * h     # implicit Euler predictive position
    x_n = copy.deepcopy(x)
    mu_lambda = BarrierEnergy.compute_mu_lambda(x, n, o, contact_area, mu)  # compute mu * lambda for each node using x^n
    # ANCHOR: dbc_initialization
    DBC_target = [] # target position of each DBC in the current time step
    for i in range(0, len(DBC)):
        if (DBC_limit[i] - x_n[DBC[i]]).dot(DBC_v[i]) > 0:
            DBC_target.append(x_n[DBC[i]] + h * DBC_v[i])
        else:
            DBC_target.append(x_n[DBC[i]])
    # ANCHOR_END: dbc_initialization

    # Newton loop
    iter = 0
    E_last = IP_val(x, e, x_tilde, m, l2, k, n, o, contact_area, (x - x_n) / h, mu_lambda, DBC, DBC_target, DBC_stiff[0], h)
    # ANCHOR: convergence_criteria
    [p, DBC_satisfied] = search_dir(x, e, x_tilde, m, l2, k, n, o, contact_area, (x - x_n) / h, mu_lambda, is_DBC, DBC, DBC_target, DBC_stiff[0], tol, h)
    while True:
        if not np.isfinite(x).all() or not np.isfinite(p).all() or not np.isfinite(E_last):
            reason = "nonFinite"
            break
        dbc_ok = sum(DBC_satisfied) == len(DBC)
        step_ok = LA.norm(p, inf) / h <= tol and dbc_ok
        raw = IP_grad(x, e, x_tilde, m, l2, k, n, o, contact_area, (x-x_n)/h, mu_lambda, DBC, DBC_target, DBC_stiff[0], h)
        free = np.asarray(raw)[~np.asarray(is_DBC, dtype=bool)]
        # 自由节点没有DBC罚梯度；动态IP梯度除h²得到离散平衡力残差。
        force = np.max(np.abs(free)) / (h*h) if free.size else 0.0
        if not np.isfinite(force):
            reason = "nonFinite"
            break
        actual_ok = small_count >= actual_step_count and dbc_ok
        met = ((step_ok if stop_mode == "step" else
                force / grad_scale <= grad_tol_rel and dbc_ok if stop_mode == "grad" else
                actual_ok if stop_mode == "actual" else
                step_ok or (actual_ok and force / grad_scale <= grad_tol_rel)))
        if met:
            reason = stop_mode
            break
        if max_iter is not None and iter >= max_iter:
            reason = "maxIter"
            break
        print('Iteration', iter, ':')
        print('residual =', LA.norm(p, inf) / h)

        if (LA.norm(p, inf) / h <= tol) & (sum(DBC_satisfied) != len(DBC)):
            # increase DBC stiffness and recompute energy value record
            DBC_stiff[0] *= 2
            E_last = IP_val(x, e, x_tilde, m, l2, k, n, o, contact_area, (x - x_n) / h, mu_lambda, DBC, DBC_target, DBC_stiff[0], h)
        # ANCHOR_END: convergence_criteria

        # filter line search
        alpha = BarrierEnergy.init_step_size(x, n, o, p)  # avoid interpenetration and tunneling
        while IP_val(x + alpha * p, e, x_tilde, m, l2, k, n, o, contact_area, (x + alpha * p - x_n) / h, mu_lambda, DBC, DBC_target, DBC_stiff[0], h) > E_last:
            alpha /= 2
        print('step size =', alpha)

        previous_x = x.copy()
        x += alpha * p
        small_count = small_count + 1 if np.max(np.abs(x-previous_x)) <= actual_step_tol else 0
        E_last = IP_val(x, e, x_tilde, m, l2, k, n, o, contact_area, (x - x_n) / h, mu_lambda, DBC, DBC_target, DBC_stiff[0], h)
        [p, DBC_satisfied] = search_dir(x, e, x_tilde, m, l2, k, n, o, contact_area, (x - x_n) / h, mu_lambda, is_DBC, DBC, DBC_target, DBC_stiff[0], tol, h)
        iter += 1

    if report is not None:
        report.update(converged=reason in ("step", "grad", "actual", "combined"), stop_reason=reason, iterations=iter, force_residual_inf_N=float(force))
    v = (x - x_n) / h   # implicit Euler velocity update
    return [x, v]

def IP_val(x, e, x_tilde, m, l2, k, n, o, contact_area, v, mu_lambda, DBC, DBC_target, DBC_stiff, h):
    return InertiaEnergy.val(x, x_tilde, m) + h * h * (     # implicit Euler
        MassSpringEnergy.val(x, e, l2, k) + 
        GravityEnergy.val(x, m) + 
        UniversalGravitationEnergy.val(x, m) + 
        BarrierEnergy.val(x, n, o, contact_area) + 
        FrictionEnergy.val(v, mu_lambda, h, n)
    ) + SpringEnergy.val(x, m, DBC, DBC_target, DBC_stiff)

def IP_grad(x, e, x_tilde, m, l2, k, n, o, contact_area, v, mu_lambda, DBC, DBC_target, DBC_stiff, h):
    return InertiaEnergy.grad(x, x_tilde, m) + h * h * (    # implicit Euler
        MassSpringEnergy.grad(x, e, l2, k) + 
        GravityEnergy.grad(x, m) + 
        UniversalGravitationEnergy.grad(x, m) + 
        BarrierEnergy.grad(x, n, o, contact_area) + 
        FrictionEnergy.grad(v, mu_lambda, h, n)
    ) + SpringEnergy.grad(x, m, DBC, DBC_target, DBC_stiff)

def IP_hess(x, e, x_tilde, m, l2, k, n, o, contact_area, v, mu_lambda, DBC, DBC_target, DBC_stiff, h):
    IJV_In = InertiaEnergy.hess(x, x_tilde, m)
    IJV_MS = MassSpringEnergy.hess(x, e, l2, k)
    IJV_B = BarrierEnergy.hess(x, n, o, contact_area)
    IJV_F = FrictionEnergy.hess(v, mu_lambda, h, n)
    IJV_S = SpringEnergy.hess(x, m, DBC, DBC_target, DBC_stiff)
    IJV_UV = UniversalGravitationEnergy.hess(x, m)
    IJV_MS[2] *= h * h    # implicit Euler
    IJV_UV[2] *= h * h    # implicit Euler
    IJV_B[2] *= h * h     # implicit Euler
    IJV_F[2] *= h * h     # implicit Euler
    IJV_In_MS = np.append(IJV_In, IJV_MS, axis=1)
    IJV_In_MS_UV = np.append(IJV_In_MS, IJV_UV, axis=1)
    IJV_In_MS_UV_B = np.append(IJV_In_MS_UV, IJV_B, axis=1)
    IJV_In_MS_UV_B_F = np.append(IJV_In_MS_UV_B, IJV_F, axis=1)
    IJV = np.append(IJV_In_MS_UV_B_F, IJV_S, axis=1)
    H = sparse.coo_matrix((IJV[2], (IJV[0], IJV[1])), shape=(len(x) * 2, len(x) * 2)).tocsr()
    return H

def search_dir(x, e, x_tilde, m, l2, k, n, o, contact_area, v, mu_lambda, is_DBC, DBC, DBC_target, DBC_stiff, tol, h):
    projected_hess = IP_hess(x, e, x_tilde, m, l2, k, n, o, contact_area, v, mu_lambda, DBC, DBC_target, DBC_stiff, h)
    reshaped_grad = IP_grad(x, e, x_tilde, m, l2, k, n, o, contact_area, v, mu_lambda, DBC, DBC_target, DBC_stiff, h).reshape(len(x) * 2, 1)
    # ANCHOR: dbc_check
    # check whether each DBC is satisfied
    DBC_satisfied = [False] * len(x)
    for i in range(0, len(DBC)):
        if LA.norm(x[DBC[i]] - DBC_target[i]) / h < tol:
            DBC_satisfied[DBC[i]] = True
    # ANCHOR_END: dbc_check
    # ANCHOR: dof_elimination
    # eliminate DOF if it's a satisfied DBC by modifying gradient and Hessian for DBC:
    for i, j in zip(*projected_hess.nonzero()):
        if (is_DBC[int(i / 2)] & DBC_satisfied[int(i / 2)]) | (is_DBC[int(j / 2)] & DBC_satisfied[int(j / 2)]): 
            projected_hess[i, j] = (i == j)
    for i in range(0, len(x)):
        if is_DBC[i] & DBC_satisfied[i]:
            reshaped_grad[i * 2] = reshaped_grad[i * 2 + 1] = 0.0
    return [spsolve(projected_hess, -reshaped_grad).reshape(len(x), 2), DBC_satisfied]
    # ANCHOR_END: dof_elimination