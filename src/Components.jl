module Components

export Plus, Minus

module Plus
using Niokr.ParametersUtils: Parameters
function A(p::Parameters, V1, V2, Vz)
    res = zeros(3, 3)

    # Первый столбец
    res[1, 1] = -p.Gamma_L - (V2 * conj(V2) * G2(p, Vz)) / p.gamma_r -
        p.Gamma_L - (conj(V1) * V1 * G1(p, Vz)) / p.gamma_r

    res[2, 1] = -V1 * conj(V2) * (G2(p, Vz) + im * F2(p, Vz) -
        G1(p, Vz) + im * F1(p, Vz)) / p.gamma_r

    res[3, 1] = -V2 * conj(V1) * (G2(p, Vz) - im * F2(p, Vz) -
        G1(p, Vz) - im * F1(p, Vz)) / p.gamma_r

    # Второй столбец
    res[1, 2] = -V2 * conj(V1) * (G2(p, Vz) + im * F2(p, Vz) -
        G1(p, Vz) + im * F1(p, Vz)) / (2.0 * p.gamma_r)

    res[2, 2] = -p.Gamma - im * (p.Omega1 - p.Omega2 - Vz * p.omega_hfs / p.c) -
        (conj(V2) * V2 * (G1(p, Vz) - im * F1(p, Vz)) +
        conj(V1) * V1 * (G2(p, Vz) + im * F2(p, Vz))) / p.gamma_r

    res[3, 2] = 0.0

    # Третий столбец
    res[1, 3] = -V1 * conj(V2) * (G2(p, Vz) - im * F2(p, Vz) -
        G1(p, Vz) - im * F1(p, Vz)) / (2.0 * p.gamma_r)

    res[2, 3] = 0.0

    res[3, 3] = -p.Gamma + im * (p.Omega1 - p.Omega2 - Vz * p.omega_hfs / p.c) -
        (conj(V2) * V2 * (G1(p, Vz) + im * F1(p, Vz)) +
        conj(V1) * V1 * (G2(p, Vz) - im * F2(p, Vz))) / p.gamma_r

    return res
end

function b(p::Parameters, V1, V2, Vz)
    res = zeros(3)
    res[1] = p.Gamma_L + conj(V1) * V1 * G1(p, Vz) / p.gamma_r
    res[2] = -(V1 * conj(V2) * (G1(p, Vz) - im * F1(p, Vz))) / p.gamma_r
    res[3] = -(V2 * conj(V1) * (G1(p, Vz) + im * F1(p, Vz))) / p.gamma_r
    return res
end

G1(p::Parameters, Vz) = p.gamma_r^2 / (p.gamma_r ^ 2 + (p.Omega1 - Vz * p.k1)^2)

G2(p::Parameters, Vz) = p.gamma_r^2 / (p.gamma_r ^ 2 + (p.Omega2 - Vz * p.k2)^2)

F1(p::Parameters, Vz) = p.gamma_r * (p.Omega1 - Vz * p.k1) /
    (p.gamma_r^2 + (p.Omega1^2 - Vz * p.k1)^2)

F2(p::Parameters, Vz) = p.gamma_r * (p.Omega2 - Vz * p.k2) /
    (p.gamma_r^2 + (p.Omega2^2 - Vz * p.k2)^2)

end

module Minus
using Niokr.ParametersUtils: Parameters
function A(p::Parameters, V1, V2, Vz)
    res = zeros(3, 3)

    # Первый столбец
    res[1, 1] = -p.Gamma_L - (V2 * conj(V2) * G2(p, Vz)) / p.gamma_r -
        p.Gamma_L - (conj(V1) * V1 * G1(p, Vz)) / p.gamma_r

    res[2, 1] = -V1 * conj(V2) * (G2(p, Vz) + im * F2(p, Vz) -
        G1(p, Vz) + im * F1(p, Vz)) / p.gamma_r

    res[3, 1] = -V2 * conj(V1) * (G2(p, Vz) - im * F2(p, Vz) -
        G1(p, Vz) - im * F1(p, Vz)) / p.gamma_r

    # Второй столбец
    res[1, 2] = -V2 * conj(V1) * (G2(p, Vz) + im * F2(p, Vz) -
        G1(p, Vz) + im * F1(p, Vz)) / (2.0 * p.gamma_r)

    res[2, 2] = -p.Gamma - im * (p.Omega1 - p.Omega2 + Vz * p.omega_hfs / p.c) -
        (conj(V2) * V2 * (G1(p, Vz) - im * F1(p, Vz)) +
        conj(V1) * V1 * (G2(p, Vz) + im * F2(p, Vz))) / p.gamma_r

    res[3, 2] = 0.0

    # Третий столбец
    res[1, 3] = -V1 * conj(V2) * (G2(p, Vz) - im * F2(p, Vz) -
        G1(p, Vz) - im * F1(p, Vz)) / (2.0 * p.gamma_r)

    res[2, 3] = 0.0

    res[3, 3] = -p.Gamma + im * (p.Omega1 - p.Omega2 + Vz * p.omega_hfs / p.c) -
        (conj(V2) * V2 * (G1(p, Vz) + im * F1(p, Vz)) +
        conj(V1) * V1 * (G2(p, Vz) - im * F2(p, Vz))) / p.gamma_r
    return res
end

function b(p::Parameters, V1, V2, Vz)
    res = zeros(3)
    res[1] = p.Gamma_L + conj(V1) * V1 * G1(p, Vz) / p.gamma_r
    res[2] = -(V1 * conj(V2) * (G1(p, Vz) - im * F1(p, Vz))) / p.gamma_r
    res[3] = -(V2 * conj(V1) * (G1(p, Vz) + im * F1(p, Vz))) / p.gamma_r
    return res
end

G1(p::Parameters, Vz) = p.gamma_r^2 / (p.gamma_r ^ 2 + (p.Omega1 + Vz * p.k1)^2)

G2(p::Parameters, Vz) = p.gamma_r^2 / (p.gamma_r ^ 2 + (p.Omega2 + Vz * p.k2)^2)

F1(p::Parameters, Vz) = p.gamma_r * (p.Omega1 + Vz * p.k1) /
    (p.gamma_r^2 + (p.Omega1^2 + Vz * p.k1)^2)

F2(p::Parameters, Vz) = p.gamma_r * (p.Omega2 + Vz * p.k2) /
    (p.gamma_r^2 + (p.Omega2^2 + Vz * p.k2)^2)

end
end