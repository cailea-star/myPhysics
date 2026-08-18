"""
@file: wigner_symbols.py
@author: cailea
@date: 2026-06-29
@brief: Reduce symbolic products of Wigner 3j symbols to Wigner 6j and 9j forms.
@ref: Varshalovich, D. A., Moskalev, A. N., and Khersonskii, V. K., Quantum Theory of Angular Momentum, 1988.
"""

from dataclasses import dataclass
from itertools import permutations

from sympy import KroneckerDelta, Mul, S, Wild, count_ops, expand, latex, powsimp, solve, sqrt, symbols
from sympy.physics.quantum.cg import CG, Wigner3j, Wigner6j, Wigner9j


def wigner3j_variants(j1_S, m1_S, j2_S, m2_S, j3_S, m3_S):
    """@math (j1 j2 j3; m1 m2 m3) = (-1)^J (j2 j1 j3; m2 m1 m3) = (-1)^J (j1 j2 j3; -m1 -m2 -m3)"""
    columns = ((j1_S, m1_S), (j2_S, m2_S), (j3_S, m3_S))
    j_sum_Expr = j1_S + j2_S + j3_S
    for permutation in permutations(range(3)):
        is_odd_I = sum(permutation[i_I] > permutation[j_I] for i_I in range(3) for j_I in range(i_I + 1, 3)) % 2
        phase_Expr = (-1) ** j_sum_Expr if is_odd_I else S.One
        permuted_columns = [columns[i_I] for i_I in permutation]
        yield phase_Expr, Wigner3j(*(value_S for j_S, m_S in permuted_columns for value_S in (j_S, m_S)))
        yield phase_Expr * (-1) ** j_sum_Expr, Wigner3j(*(value_S for j_S, m_S in permuted_columns for value_S in (j_S, -m_S)))


def simplify_wigner9j_phase(expr_Expr):
    """@math {9j}_odd = (-1)^Σj {9j}"""
    wigner9j_Expr = next(iter(expr_Expr.atoms(Wigner9j)), None)
    if wigner9j_Expr is None:
        return expr_Expr

    zero_conditions, integer_conditions = wigner3j_constraints(expr_Expr.atoms(Wigner3j))
    integer_conditions.append(sum(wigner9j_Expr.args))

    matrix = [wigner9j_Expr.args[i_I:i_I + 3] for i_I in range(0, 9, 3)]
    odd_variant_Expr = Wigner9j(*(matrix[i_I][j_I] for i_I in (0, 2, 1) for j_I in (2, 0, 1)))
    odd_Expr = expr_Expr.xreplace({wigner9j_Expr: (-1) ** sum(wigner9j_Expr.args) * odd_variant_Expr})
    candidates = (simplify_phase(expr_Expr, zero_conditions, integer_conditions), simplify_phase(odd_Expr, zero_conditions, integer_conditions))
    return min(candidates, key=lambda candidate_Expr: (count_ops(candidate_Expr), str(candidate_Expr)))


def wigner3j_constraints(wigner3js):
    """@math m1 + m2 + m3 = 0, j1 + j2 + j3 ∈ Z"""
    zero_conditions = []
    integer_conditions = []
    for wigner3j_Expr in wigner3js:
        j1_S, m1_S, j2_S, m2_S, j3_S, m3_S = wigner3j_Expr.args
        zero_conditions.append(m1_S + m2_S + m3_S)
        integer_conditions.extend([
            j1_S + j2_S + j3_S, j1_S + j2_S - j3_S,
            j1_S - j2_S + j3_S, j1_S - j2_S - j3_S,
            j1_S + m1_S, j1_S - m1_S,
            j2_S + m2_S, j2_S - m2_S,
            j3_S + m3_S, j3_S - m3_S,
            2 * j1_S, 2 * m1_S,
            2 * j2_S, 2 * m2_S,
            2 * j3_S, 2 * m3_S,
        ])
    return zero_conditions, integer_conditions


def _eliminate_sum_symbols(exponent_Expr, zero_conditions, sum_symbols, zero_symbols):
    """@math m1 + m2 + m3 = 0 → m_k = -Σ_{i≠k}m_i"""
    remaining_symbols = exponent_Expr.free_symbols & set(sum_symbols)
    while remaining_symbols:
        candidates = []
        for condition_Expr in zero_conditions:
            for symbol_S in condition_Expr.free_symbols & remaining_symbols:
                solutions = solve(condition_Expr, symbol_S)
                if solutions:
                    candidate_Expr = expand(exponent_Expr.subs(symbol_S, solutions[0]))
                    if len(candidate_Expr.free_symbols & set(sum_symbols)) < len(remaining_symbols): candidates.append(candidate_Expr)
        if not candidates: break
        exponent_Expr = min(candidates, key=lambda candidate_Expr: (len(candidate_Expr.free_symbols & zero_symbols), count_ops(candidate_Expr), str(candidate_Expr)))
        remaining_symbols = exponent_Expr.free_symbols & set(sum_symbols)
    return exponent_Expr


def _reduce_phase_exponent(exponent_Expr, zero_conditions, integer_conditions, zero_symbols):
    """@math (-1)^(e + 2n) = (-1)^e"""
    exponent_Expr = expand(exponent_Expr)
    solve_symbols = sorted(exponent_Expr.free_symbols & zero_symbols, key=str)
    solutions = solve(zero_conditions, solve_symbols, dict=True) if solve_symbols else []
    substitutions = solutions[0] if solutions else {}
    exponent_Expr = expand(exponent_Expr.subs(substitutions))
    integer_conditions = [expand(condition_Expr.subs(substitutions)) for condition_Expr in integer_conditions]
    symbols_S1D = sorted(set().union(exponent_Expr.free_symbols, *(condition_Expr.free_symbols for condition_Expr in integer_conditions)), key=str)
    coefficients_I1D = [int(exponent_Expr.coeff(symbol_S)) for symbol_S in symbols_S1D]
    rows_I2D = list(dict.fromkeys(tuple(int(condition_Expr.coeff(symbol_S)) % 2 for symbol_S in symbols_S1D) for condition_Expr in integer_conditions))
    rows_I2D = [list(row_I1D) for row_I1D in rows_I2D if any(row_I1D)]
    # c = (c mod 2) + 2b (mod 4), b in GF(2).
    bits_I1D = [(coefficient_I % 4) // 2 for coefficient_I in coefficients_I1D]
    pivot_I = 0
    for symbol_I in range(len(symbols_S1D)):
        row_I = next((row_I for row_I in range(pivot_I, len(rows_I2D)) if rows_I2D[row_I][symbol_I]), None)
        if row_I is None: continue
        rows_I2D[pivot_I], rows_I2D[row_I] = rows_I2D[row_I], rows_I2D[pivot_I]
        for row_I in range(len(rows_I2D)):
            if row_I != pivot_I and rows_I2D[row_I][symbol_I]: rows_I2D[row_I] = [value_I ^ pivot_value_I for value_I, pivot_value_I in zip(rows_I2D[row_I], rows_I2D[pivot_I])]
        if bits_I1D[symbol_I]: bits_I1D = [value_I ^ pivot_value_I for value_I, pivot_value_I in zip(bits_I1D, rows_I2D[pivot_I])]
        pivot_I += 1
        if pivot_I == len(rows_I2D): break
    residues_I1D = [(coefficient_I % 2 + 2 * bit_I) % 4 for coefficient_I, bit_I in zip(coefficients_I1D, bits_I1D)]
    residues_I1D = [-1 if residue_I == 3 else residue_I for residue_I in residues_I1D]
    constant_Expr = exponent_Expr.subs({symbol_S: 0 for symbol_S in symbols_S1D})
    return expand(int(constant_Expr) % 2 + sum(residue_I * symbol_S for residue_I, symbol_S in zip(residues_I1D, symbols_S1D)))


def simplify_phase(expr_Expr, zero_conditions, integer_conditions, eliminate_symbols=()):
    """@math (-1)^(e + 2n) = (-1)^e"""
    expr_Expr = powsimp(expr_Expr, force=True)
    zero_symbols = set().union(*(condition_Expr.free_symbols for condition_Expr in zero_conditions))

    def simplify_power(power_Expr):
        """@math (-1)^e = (-1)^e'"""
        exponent_Expr = _reduce_phase_exponent(power_Expr.exp, zero_conditions, integer_conditions, zero_symbols)
        exponent_Expr = _eliminate_sum_symbols(exponent_Expr, zero_conditions, eliminate_symbols, zero_symbols)
        return (-1) ** exponent_Expr

    expr_Expr = expr_Expr.replace(lambda value_Expr: value_Expr.is_Pow and value_Expr.base == -1, simplify_power)
    return powsimp(expr_Expr, force=True)


def _split_wigner3j_product(expr_Expr):
    """@math E = P Π_i Wigner3j_i"""
    wigner3js = tuple(factor_Expr for factor_Expr in expr_Expr.as_ordered_factors() if isinstance(factor_Expr, Wigner3j))
    return expr_Expr / Mul(*wigner3js), wigner3js


def _cg_to_w3j(cg_Expr):
    """@math CG(j1m1,j2m2|j3m3) = (-1)^(j1-j2+m3)√(2j3+1) Wigner3j(j1,m1,j2,m2,j3,-m3)"""
    j1_S, m1_S, j2_S, m2_S, j3_S, m3_S = cg_Expr.args
    return (-1) ** (j1_S - j2_S + m3_S) * sqrt(2 * j3_S + 1) * Wigner3j(j1_S, m1_S, j2_S, m2_S, j3_S, -m3_S)


@dataclass(frozen=True)
class WignerReductionPattern:
    sum_symbols: tuple
    lhs_Expr: object
    rhs_Expr: object

    def match_wigner3j_product(self, expr_w3js, pattern_w3js, matched_symbols=None, phase_Expr=S.One):
        """@math Π_i W_i = (-1)^e Π_i W'_i"""
        matched_symbols = matched_symbols or {}
        if not pattern_w3js:
            yield phase_Expr, matched_symbols, expr_w3js
            return

        pattern_w3j_Expr = pattern_w3js[0]
        rest_pattern_w3js = pattern_w3js[1:]
        for w3j_I, expr_w3j_Expr in enumerate(expr_w3js):
            rest_expr_w3js = expr_w3js[:w3j_I] + expr_w3js[w3j_I + 1:]
            for variant_phase_Expr, variant_Expr in wigner3j_variants(*pattern_w3j_Expr.args):
                match = expr_w3j_Expr.match(variant_Expr)
                if match is None:
                    continue
                if any(matched_symbols.get(key_S, value_Expr) != value_Expr for key_S, value_Expr in match.items()):
                    continue
                yield from self.match_wigner3j_product(rest_expr_w3js, rest_pattern_w3js, matched_symbols | match, phase_Expr * variant_phase_Expr)

    def apply_wigner_reduction(self, sum_symbols_in, expr_in_Expr):
        """@math Σ_dummy lhs = rhs"""
        expr_w3j_Expr = expr_in_Expr.replace(lambda value_Expr: isinstance(value_Expr, CG), _cg_to_w3j)
        wild_symbols = {symbol_S: Wild(str(symbol_S)) for symbol_S in self.lhs_Expr.free_symbols | self.rhs_Expr.free_symbols}
        lhs_wild_Expr = self.lhs_Expr.xreplace(wild_symbols)
        rhs_wild_Expr = self.rhs_Expr.xreplace(wild_symbols)

        expr_rest_Expr, expr_w3js = _split_wigner3j_product(expr_w3j_Expr)
        lhs_rest_Expr, lhs_w3js = _split_wigner3j_product(lhs_wild_Expr)
        sum_wilds = tuple(wild_symbols[symbol_S] for symbol_S in self.sum_symbols)
        sum_symbol_set = set(sum_symbols_in)

        # E = E_rest E_matched.
        if len(expr_w3js) < len(lhs_w3js):
            return sum_symbols_in, expr_in_Expr

        zero_conditions, integer_conditions = wigner3j_constraints(expr_w3js)
        for matched_phase_Expr, matched_symbols, unmatched_w3js in self.match_wigner3j_product(expr_w3js, lhs_w3js):
            reduced_symbols = {-matched_symbols[wild_S] if matched_symbols[wild_S].could_extract_minus_sign() else matched_symbols[wild_S] for wild_S in sum_wilds}
            if len(reduced_symbols) != len(sum_wilds) or not reduced_symbols <= sum_symbol_set:
                continue

            prefactor_Expr = expr_rest_Expr * Mul(*unmatched_w3js) / (lhs_rest_Expr.xreplace(matched_symbols) * matched_phase_Expr.xreplace(matched_symbols))
            reduced_Expr = powsimp(prefactor_Expr * rhs_wild_Expr.xreplace(matched_symbols), force=True)
            reduced_Expr = simplify_phase(reduced_Expr, zero_conditions, integer_conditions, reduced_symbols)
            if reduced_Expr.free_symbols & reduced_symbols:
                continue
            return tuple(symbol_S for symbol_S in sum_symbols_in if symbol_S not in reduced_symbols), reduced_Expr

        return sum_symbols_in, expr_in_Expr


a_S, alpha_S = symbols("a alpha")
ap_S, alphap_S = symbols("ap alphap")
b_S, beta_S = symbols("b beta")
c_S, gamma_S = symbols("c gamma")
p_S, psi_S = symbols("p psi")
q_S, kappa_S = symbols("q kappa")
r_S, rho_S = symbols("r rho")
s_S, sigma_S = symbols("s sigma")
t_S, tau_S = symbols("t tau")
u_S, nu_S = symbols("u nu")


THREE_3J_SUM_3_TO_ONE_6J_ONE_3J = WignerReductionPattern(
    sum_symbols=(kappa_S, psi_S, rho_S),
    lhs_Expr=(
        (-1) ** (p_S - psi_S + q_S - kappa_S + r_S - rho_S)
        * Wigner3j(p_S, psi_S, a_S, alpha_S, q_S, -kappa_S)
        * Wigner3j(q_S, kappa_S, b_S, beta_S, r_S, -rho_S)
        * Wigner3j(r_S, rho_S, c_S, gamma_S, p_S, -psi_S)
    ),
    rhs_Expr=(
        1 
        * Wigner3j(a_S, -alpha_S, b_S, -beta_S, c_S, -gamma_S)
        * Wigner6j(a_S, b_S, c_S, r_S, p_S, q_S)
    ),
)

FOUR_3J_SUM_5_TO_ONE_6J = WignerReductionPattern(
    sum_symbols=(kappa_S, psi_S, rho_S, sigma_S, tau_S),
    lhs_Expr=(
        (-1) ** (p_S - psi_S + q_S - kappa_S + r_S - rho_S + s_S - sigma_S + t_S - tau_S)
        * Wigner3j(p_S, psi_S, a_S, -alpha_S, q_S, kappa_S)
        * Wigner3j(q_S, -kappa_S, r_S, rho_S, t_S, tau_S)
        * Wigner3j(r_S, -rho_S, ap_S, alphap_S, s_S, sigma_S)
        * Wigner3j(s_S, -sigma_S, p_S, -psi_S, t_S, -tau_S)
    ),
    rhs_Expr=(
        (-1) ** (a_S - alpha_S) / (2 * a_S + 1)
        * Wigner6j(q_S, p_S, a_S, s_S, r_S, t_S)
        * KroneckerDelta(a_S, ap_S)
        * KroneckerDelta(alpha_S, alphap_S)
    ),
)


FIVE_3J_SUM_6_TO_ONE_9J_ONE_3J = WignerReductionPattern(
    sum_symbols=(psi_S, kappa_S, rho_S, tau_S, sigma_S, nu_S),
    lhs_Expr=(
        (-1) ** (p_S - psi_S + q_S - kappa_S + r_S - rho_S + s_S - sigma_S + t_S - tau_S + u_S - nu_S)
        * Wigner3j(p_S, psi_S, a_S, alpha_S, q_S, kappa_S)
        * Wigner3j(q_S, -kappa_S, t_S, -tau_S, r_S, -rho_S)
        * Wigner3j(r_S, rho_S, b_S, beta_S, s_S, sigma_S)
        * Wigner3j(s_S, -sigma_S, p_S, -psi_S, u_S, -nu_S)
        * Wigner3j(u_S, nu_S, c_S, gamma_S, t_S, tau_S)
    ),
    rhs_Expr=(
        (-1) ** (r_S + b_S + s_S)
        * Wigner3j(a_S, -alpha_S, b_S, -beta_S, c_S, -gamma_S)
        * Wigner9j(a_S, b_S, c_S, p_S, s_S, u_S, q_S, r_S, t_S)
    ),
)

def _print_latex_case(case_text, input_text, reference_Expr, computed_Expr, latex_symbol_names):
    """@math input -> reference = computed"""
    print()
    for label, body in (("Case", case_text), ("Input LaTeX", input_text), ("Reference LaTeX", latex(reference_Expr, symbol_names=latex_symbol_names)), ("Computed LaTeX", latex(computed_Expr, symbol_names=latex_symbol_names))):
        print(f"[{label}]\n$$\n{body}\n$$")

if __name__ == "__main__":
    s_S, ms_S, sp_S, msp_S = symbols("s ms sp msp")
    l_S, ml_S, lp_S, mlp_S = symbols("l ml lp mlp")
    j_S, mj_S, jp_S, mjp_S = symbols("j mj jp mjp")
    lamb_S, mu_S = symbols("lamb mu")
    lambs_S, mus_S, lambl_S, mul_S = symbols("lambs mus lambl mul")
    Bred_S = symbols("Bred")
    latex_symbol_names = {
        ms_S: r"m_s",
        sp_S: r"s'",
        msp_S: r"m_{s'}",
        ml_S: r"m_l",
        lp_S: r"l'",
        mlp_S: r"m_{l'}",
        mj_S: r"m_j",
        jp_S: r"j'",
        mjp_S: r"m_{j'}",
        lamb_S: r"\lambda",
        mu_S: r"\mu",
        lambs_S: r"\lambda_s",
        mus_S: r"\mu_s",
        lambl_S: r"\lambda_l",
        mul_S: r"\mu_l",
        Bred_S: r"\langle l\|B^{(\lambda)}\|l'\rangle",
    }

    scalar_input_Expr = (
        CG(s_S, ms_S, l_S, ml_S, j_S, mj_S)
        * CG(sp_S, msp_S, lp_S, mlp_S, jp_S, mjp_S)
        * (-1) ** (l_S - ml_S)
        * Wigner3j(l_S, -ml_S, lamb_S, -mu_S, lp_S, mlp_S)
        * (-1) ** (s_S - ms_S)
        * Wigner3j(s_S, -ms_S, lamb_S, mu_S, sp_S, msp_S)
        * (-1) ** (lamb_S - mu_S)
        / sqrt(2 * lamb_S + 1)
    )
    scalar_sum_symbols = (ms_S, ml_S, msp_S, mlp_S, mu_S)
    scalar_reduced_symbols, scalar_reduced_Expr = FOUR_3J_SUM_5_TO_ONE_6J.apply_wigner_reduction(scalar_sum_symbols, scalar_input_Expr)
    scalar_reference_Expr = (-1) ** (jp_S + l_S + lamb_S + sp_S) * sqrt(2 * jp_S + 1) * KroneckerDelta(j_S, jp_S) * KroneckerDelta(mj_S, mjp_S) * Wigner6j(l_S, s_S, j_S, sp_S, lp_S, lamb_S) / sqrt((2 * j_S + 1) * (2 * lamb_S + 1))
    _print_latex_case(r"\langle (sl)jm_j | [A_\lambda(s)B_\lambda(l)]_{00} | (s'l')j'm_j' \rangle", r"\sum_{" + ",".join(latex(symbol_S, symbol_names=latex_symbol_names) for symbol_S in scalar_sum_symbols) + r"} \langle j m_j\mid s m_s,l m_l\rangle\langle s'm_{s'},l'm_{l'}\mid j'm_{j'}\rangle\;" + latex(scalar_input_Expr / (CG(s_S, ms_S, l_S, ml_S, j_S, mj_S) * CG(sp_S, msp_S, lp_S, mlp_S, jp_S, mjp_S)), symbol_names=latex_symbol_names), scalar_reference_Expr, scalar_reduced_Expr, latex_symbol_names)
    assert scalar_reduced_symbols == () and scalar_reduced_Expr == scalar_reference_Expr

    tensor_input_Expr = (
        CG(s_S, ms_S, l_S, ml_S, j_S, mj_S)
        * CG(sp_S, msp_S, lp_S, mlp_S, jp_S, mjp_S)
        * CG(lambs_S, mus_S, lambl_S, mul_S, lamb_S, mu_S)
        * (-1) ** (l_S - ml_S)
        * Wigner3j(l_S, -ml_S, lambl_S, mul_S, lp_S, mlp_S)
        * (-1) ** (s_S - ms_S)
        * Wigner3j(s_S, -ms_S, lambs_S, mus_S, sp_S, msp_S)
    )
    tensor_sum_symbols = (ms_S, ml_S, msp_S, mlp_S, mus_S, mul_S)
    tensor_reduced_symbols, tensor_reduced_Expr = FIVE_3J_SUM_6_TO_ONE_9J_ONE_3J.apply_wigner_reduction(tensor_sum_symbols, tensor_input_Expr)
    tensor_reduced_Expr = simplify_wigner9j_phase(tensor_reduced_Expr)
    tensor_reference_Expr = (-1) ** (j_S - mj_S) * sqrt((2 * j_S + 1) * (2 * jp_S + 1) * (2 * lamb_S + 1)) * Wigner3j(lamb_S, mu_S, jp_S, mjp_S, j_S, -mj_S) * Wigner9j(j_S, lamb_S, jp_S, l_S, lambl_S, lp_S, s_S, lambs_S, sp_S)
    _print_latex_case(r"\langle (sl)jm_j | [A_{\lambda_s}(s)B_{\lambda_l}(l)]_{\lambda\mu} | (s'l')j'm_j' \rangle", r"\sum_{" + ",".join(latex(symbol_S, symbol_names=latex_symbol_names) for symbol_S in tensor_sum_symbols) + r"} \langle j m_j\mid s m_s,l m_l\rangle\langle s'm_{s'},l'm_{l'}\mid j'm_{j'}\rangle\langle\lambda_s\mu_s,\lambda_l\mu_l\mid\lambda\mu\rangle\;" + latex(tensor_input_Expr / (CG(s_S, ms_S, l_S, ml_S, j_S, mj_S) * CG(sp_S, msp_S, lp_S, mlp_S, jp_S, mjp_S) * CG(lambs_S, mus_S, lambl_S, mul_S, lamb_S, mu_S)), symbol_names=latex_symbol_names), tensor_reference_Expr, tensor_reduced_Expr, latex_symbol_names)
    assert tensor_reduced_symbols == () and tensor_reduced_Expr == tensor_reference_Expr

    subsystem_input_Expr = (
        CG(s_S, ms_S, l_S, ml_S, j_S, mj_S)
        * CG(s_S, ms_S, lp_S, mlp_S, jp_S, mjp_S)
        * (-1) ** (l_S - ml_S)
        * Wigner3j(l_S, -ml_S, lamb_S, mu_S, lp_S, mlp_S)
        * KroneckerDelta(s_S, sp_S)
        * Bred_S
    )
    subsystem_sum_symbols = (ms_S, ml_S, mlp_S)
    subsystem_reduced_symbols, subsystem_reduced_Expr = THREE_3J_SUM_3_TO_ONE_6J_ONE_3J.apply_wigner_reduction(subsystem_sum_symbols, subsystem_input_Expr)
    subsystem_reduced_Expr = subsystem_reduced_Expr.xreplace({Wigner3j(lamb_S, -mu_S, jp_S, -mjp_S, j_S, mj_S): (-1) ** (lamb_S + jp_S + j_S) * Wigner3j(j_S, -mj_S, lamb_S, mu_S, jp_S, mjp_S), Wigner6j(lamb_S, jp_S, j_S, s_S, l_S, lp_S): Wigner6j(l_S, j_S, s_S, jp_S, lp_S, lamb_S)})
    subsystem_reduced_Expr = simplify_phase(powsimp(subsystem_reduced_Expr, force=True), *wigner3j_constraints(subsystem_reduced_Expr.atoms(Wigner3j)))
    subsystem_reference_Expr = (
        (-1) ** (j_S - mj_S)
        * Wigner3j(j_S, -mj_S, lamb_S, mu_S, jp_S, mjp_S)
        * KroneckerDelta(s_S, sp_S)
        * (-1) ** (s_S + lp_S + j_S + lamb_S)
        * sqrt((2 * j_S + 1) * (2 * jp_S + 1))
        * Wigner6j(l_S, j_S, s_S, jp_S, lp_S, lamb_S)
        * Bred_S
    )
    _print_latex_case(r"\langle (sl)jm_j | B_{\lambda\mu}(l) | (s'l')j'm_j' \rangle", r"\sum_{" + ",".join(latex(symbol_S, symbol_names=latex_symbol_names) for symbol_S in subsystem_sum_symbols) + r"} \langle j m_j\mid s m_s,l m_l\rangle\langle s m_s,l'm_{l'}\mid j'm_{j'}\rangle\;" + latex(subsystem_input_Expr / (CG(s_S, ms_S, l_S, ml_S, j_S, mj_S) * CG(s_S, ms_S, lp_S, mlp_S, jp_S, mjp_S)), symbol_names=latex_symbol_names), subsystem_reference_Expr, subsystem_reduced_Expr, latex_symbol_names)
    assert subsystem_reduced_symbols == () and simplify_phase(powsimp(subsystem_reduced_Expr / subsystem_reference_Expr, force=True), *wigner3j_constraints(subsystem_reference_Expr.atoms(Wigner3j))) == 1
