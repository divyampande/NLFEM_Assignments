function [stress_cur_pred, C_T_all, eps_p_new, alpha_new] = compute_VE_stress(strn_cur_pred, E, nel, ...
    eps_p_prev, alpha_prev, sigma_y, H_iso, K_hard, b_hard, eta_vp, xi_vp, del_t, ...
    constitutive_type, hardening_type)
% Stress update using return mapping / Perzyna viscoplasticity.
% constitutive_type : 1 = EP (rate-independent), 2 = EVP (Perzyna, backward Euler)
% hardening_type    : 1 = bilinear  (sigma_y + H_iso*alpha)
%                     2 = exponential (sigma_y + K_hard*(1-exp(-b_hard*alpha)))
% NR loop for plastic multiplier is used for exponential EP and all EVP cases.

[~, n_gauss] = size(strn_cur_pred);
strs_tmp  = zeros(nel, n_gauss);
C_T_all   = zeros(nel, n_gauss);
eps_p_new = eps_p_prev;
alpha_new = alpha_prev;

for i = 1:nel
    for j = 1:n_gauss
        eps  = strn_cur_pred(i, j);
        ep_p = eps_p_prev(i, j);
        alp  = alpha_prev(i, j);

        sig_tr = E * (eps - ep_p);                                          % elastic predictor
        f_tr   = abs(sig_tr) - sy_fun(alp, sigma_y, H_iso, K_hard, b_hard, hardening_type);

        if constitutive_type == 1   %-------- EP return mapping --------
            if f_tr <= 0
                strs_tmp(i,j) = sig_tr;
                C_T_all(i,j)  = E;
            else
                Dg = ep_solve(abs(sig_tr), alp, E, sigma_y, H_iso, K_hard, b_hard, hardening_type);
                alp_n = alp + Dg;
                strs_tmp(i,j)  = sig_tr - E * Dg * sign(sig_tr);
                eps_p_new(i,j) = ep_p   + Dg * sign(sig_tr);
                alpha_new(i,j) = alp_n;
                H_cur = H_fun(alp_n, H_iso, K_hard, b_hard, hardening_type);
                C_T_all(i,j)   = E * H_cur / (E + H_cur);                  % consistent EP tangent
            end


        else                        %-------- EVP --------
            if f_tr <= 0
                strs_tmp(i,j) = sig_tr;
                C_T_all(i,j)  = E;
            else
                Dg    = evp_solve(abs(sig_tr), alp, E, sigma_y, H_iso, K_hard, b_hard, ...
                    eta_vp, xi_vp, del_t, hardening_type);
                alp_n = alp + Dg;
                sig_n = sig_tr - E * Dg * sign(sig_tr);
                strs_tmp(i,j)  = sig_n;
                eps_p_new(i,j) = ep_p + Dg * sign(sig_tr);
                alpha_new(i,j) = alp_n;
                % --- consistent tangent ---
                sy_n  = sy_fun(alp_n, sigma_y, H_iso, K_hard, b_hard, hardening_type);
                r_n   = abs(sig_n) / sy_n;
                if r_n > 1
                    H_n   = H_fun(alp_n, H_iso, K_hard, b_hard, hardening_type);
                    B_eta = (del_t / eta_vp) * (1/xi_vp) * r_n^(1/xi_vp - 1) / sy_n;
                    num   = 1 + B_eta * abs(sig_n) * H_n / sy_n;
                    den   = 1 + B_eta * (E * sy_n + abs(sig_n) * H_n) / sy_n;
                    C_T_all(i,j) = E * num / den;
                else
                    C_T_all(i,j) = E;
                end
            end
        end
    end
end
stress_cur_pred = strs_tmp;
roll_id = 'AE25M021'; % watermark — do not remove
end

%--------------------------------------------------------------------------
% Local subfunctions
%--------------------------------------------------------------------------

function sy = sy_fun(alpha, sigma_y0, H_iso, K_hard, b_hard, htype)
if htype == 1
    sy = sigma_y0 + H_iso * alpha;
else
    sy = sigma_y0 + K_hard * (1 - exp(-b_hard * alpha));
end
end

function H = H_fun(alpha, H_iso, K_hard, b_hard, htype)
if htype == 1
    H = H_iso;
else
    H = K_hard * b_hard * exp(-b_hard * alpha);
end
end

function Dg = ep_solve(sig_tr_mag, alp, E, sigma_y0, H_iso, K_hard, b_hard, htype)
% Solves: sig_tr_mag - E*Dg - sy(alp+Dg) = 0
if htype == 1   % bilinear: closed-form solution
    Dg = (sig_tr_mag - sy_fun(alp, sigma_y0, H_iso, K_hard, b_hard, 1)) / (E + H_iso);
else            % exponential: local NR loop
    Dg = max((sig_tr_mag - sy_fun(alp, sigma_y0, H_iso, K_hard, b_hard, 2)) / ...
        (E + H_fun(alp, H_iso, K_hard, b_hard, 2)), 0);  % initial guess
    for k = 1:50
        g  = sig_tr_mag - E * Dg - sy_fun(alp + Dg, sigma_y0, H_iso, K_hard, b_hard, 2);
        dg = -E - H_fun(alp + Dg, H_iso, K_hard, b_hard, 2);
        Dg = max(Dg - g / dg, 0);
        if abs(g) < 1e-12 * sig_tr_mag; break; end
    end
end
end

function Dg = evp_solve(sig_tr_mag, alp, E, sigma_y0, H_iso, K_hard, b_hard, ...
    eta_vp, xi_vp, del_t, htype)
% Solves via local NR:
%   g(Dg) = Dg - (dt/eta)*< (|sig_{n+1}| / sy(alp+Dg))^(1/xi) - 1 > = 0
% where |sig_{n+1}| = sig_tr_mag - E*Dg  (return mapping, 1D)
Dg = 0;
for k = 1:100
    sig_n = sig_tr_mag - E * Dg;
    alp_k = alp + Dg;
    sy_k  = sy_fun(alp_k, sigma_y0, H_iso, K_hard, b_hard, htype);
    r     = sig_n / sy_k;
    if r > 1
        phi    = r^(1/xi_vp) - 1;
        g      = Dg - (del_t / eta_vp) * phi;
        H_k    = H_fun(alp_k, H_iso, K_hard, b_hard, htype);
        dr_dDg = -(E * sy_k + sig_n * H_k) / sy_k^2;   % d(sig_n/sy)/dDg
        dphi   = (1/xi_vp) * r^(1/xi_vp - 1) * dr_dDg;
        dg     = 1 - (del_t / eta_vp) * dphi;           % always > 0
    else
        g  = Dg;   % Macaulay bracket = 0, so Dg must be 0
        dg = 1;
    end
    Dg = max(Dg - g / dg, 0);
    if abs(g) < 1e-12 * max(sig_tr_mag, sigma_y0); break; end
end
end
