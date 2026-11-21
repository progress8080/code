"""Unit conversions. Every calculation inside the package is in SI (N, m, Pa)."""

G = 9.80665                      # m/s², standard gravity

TONNE_FORCE = 1000.0 * G         # N per tonne (mass) under standard gravity
MM = 1e-3                        # m per mm
MPA = 1e6                        # Pa per MPa
GPA = 1e9                        # Pa per GPa
KN_PER_MM = 1e3 / MM             # N/m per kN/mm


def to_tonnes(force_n):
    """Force in N -> tonnes (also N·m -> t·m, N/m -> t/m)."""
    return force_n / TONNE_FORCE
