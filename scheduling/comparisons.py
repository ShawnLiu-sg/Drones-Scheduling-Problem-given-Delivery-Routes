"""Arithmetic on certified paired delay bounds."""
PAIR_STATES = ("严格变化", "零变化", "未定", "至少一端已证不可行")
def subtract(a, b):
    return a - b if a is not None and b is not None else None


def pair_status(a, b):
    """The more restrictive setting is A. Keep raw and monotone bounds separately."""
    la, ua, lb, ub = a["lower"], a["upper"], b["lower"], b["upper"]
    sa, sb = a["status"], b["status"]
    raw_low, raw_high = subtract(la, ub), subtract(ua, lb)
    if sb == "infeasible":
        assert ua is None
        return raw_low, raw_high, None, None, None, PAIR_STATES[3], "两端不可行（可行域嵌套）"
    if sa == "infeasible":
        evidence = "严格端不可行，宽松端可行" if ub is not None else "严格端不可行，宽松端未定"
        return raw_low, raw_high, None, None, None, PAIR_STATES[3], evidence
    low = max(0, raw_low or 0)
    high = raw_high
    assert high is None or low <= high
    # A finite witness for the restrictive setting also certifies feasibility
    # of the relaxation. A lower bound alone leaves feasibility at A unknown.
    if ua is None:
        evidence = "严格端若可行则严格变化" if low > 0 else "严格端可行性未定"
        return raw_low, raw_high, low, high, None, PAIR_STATES[2], evidence
    exact = low if high == low else None
    if low > 0:
        state = PAIR_STATES[0]
    elif high == 0:
        state = PAIR_STATES[1]
    else:
        state = PAIR_STATES[2]
    evidence = "两端已有可行解" if ub is not None else "宽松端由可行域嵌套证实可行"
    return raw_low, raw_high, low, high, exact, state, evidence


def triple_bounds(bounds):
    """Sharp linear contrast bounds given endpoint intervals and Z1 >= Z2 >= Z3."""
    e1, e2, e3 = bounds
    if any(e["status"] == "infeasible" for e in (e1, e2, e3)):
        return None, None, None, None, "无法判断", "至少一端已证不可行"
    l1, l2, l3 = (e["lower"] for e in bounds)
    u1, u2, u3 = (e["upper"] for e in bounds)
    raw_low = None if u2 is None else l1 - 2 * u2 + l3
    raw_high = None if u1 is None or u3 is None else u1 - 2 * l2 + u3
    finite_z2_max = [v for v in (u1, u2) if v is not None]
    z2_max = min(finite_z2_max) if finite_z2_max else None
    z2_min = max(l2, l3)
    low = None if z2_max is None else max(l1, z2_max) - 2 * z2_max + l3
    z3_max = min(u3, z2_min) if u3 is not None else z2_min
    high = None if u1 is None else u1 - 2 * z2_min + z3_max
    if low is not None and high is not None:
        assert low <= high
    if u1 is None:
        return raw_low, raw_high, low, high, "无法判断", "最紧端无可行解记录"
    if low >= 0:
        strength = "严格递减" if low > 0 else ("相等" if high == 0 else "可证非增")
        return raw_low, raw_high, low, high, "支持非增", strength
    if high < 0:
        return raw_low, raw_high, low, high, "不支持", "第三架严格优于第二架的边际收益"
    return raw_low, raw_high, low, high, "无法判断", "差值区间不能判定"
