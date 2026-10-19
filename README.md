![MORDOR logo](assets/MORDOR_LOGO.png)

# MORDOR

**MOdel for Redox-Dependent Outgassing Research**

MORDOR is a Python model for calculating volatile degassing and gas speciation in the **C–H–O–N–S system**. It combines chemical equilibrium, volatile mass balances, and melt solubility laws to investigate how pressure, temperature, oxygen fugacity, melt composition, and volatile inventory affect the gases released from silicate melts.

The model includes separate formulations for intrusive and extrusive degassing, optional volatile partitioning during partial melting, and a choice of solubility parameterizations.

## Model capabilities

- Calculates ten gas species: **H₂, H₂O, CH₄, CO, CO₂, S₂, SO₂, H₂S, N₂, and NH₃**.
- Returns dissolved H₂O, CO₂, sulfur, and nitrogen abundances in the melt.
- Defines oxygen fugacity relative to the **IW**, **QFM**, or **NiNiO** buffer, with a user-specified logarithmic offset.
- Offers ideal-gas calculations and approximate fugacity corrections for selected species.
- Allows direct specification of melt volatile contents or calculation of rock-to-melt partitioning during batch or fractional melting.
- Includes an Fe–FeO–Fe₂O₃ equilibrium calculation to account for the redox-dependent distribution of iron.

Each solver call evaluates a single set of conditions. Pressure paths, parameter sweeps, and atmospheric evolution must be constructed by calling the solvers from an external script.

## Installation

MORDOR requires **Python 3**, **NumPy**, **SciPy**, and **Matplotlib**. Matplotlib is imported by supporting modules, including when no plots are produced.

Clone the repository and enter its directory:

```bash
git clone https://github.com/CBrachmann/MORDOR---MOdel-for-Redox-Dependent-Outgassing-Research.git
cd MORDOR---MOdel-for-Redox-Dependent-Outgassing-Research
```

Optionally create a virtual environment:

```bash
python -m venv .venv
```

Activate it on Linux or macOS:

```bash
source .venv/bin/activate
```

Or in Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the dependencies:

```bash
python -m pip install numpy scipy matplotlib
```

The repository is organized as Python scripts rather than an installable package. Run examples from the repository directory so that the local modules can be imported.

## Quick start

1. Open [`Set_up.py`](Set_up.py).
2. Set the melt composition, volatile inventory, pressure, temperature, redox buffer, and solubility laws.
3. Choose the degassing style and whether to calculate volatile partitioning.
4. Run the script:

```bash
python Set_up.py
```

`Set_up.py` stores the results in Python variables. It does **not** print results, generate plots, or save output files automatically. To inspect the results interactively, run:

```bash
python -i Set_up.py
```

For the default intrusive calculation, inspect variables such as:

```python
print("CO2 gas mole fraction:", XCO2)
print("Dissolved H2O mole fraction:", XH2O_melt)
print("Gas-phase parameter alphaG:", alphaG)
```

The supplied setup uses a basaltic composition, a temperature of **1675 K**, a pressure of **100 bar**, an oxygen fugacity of **IW + 3**, ideal-gas behavior, and directly specified melt volatile contents.

## Inputs and units

### General settings

| Variable in `Set_up.py` | Meaning | Values or units |
| --- | --- | --- |
| `T` | Melt temperature | Kelvin |
| `P` | Melt pressure for intrusive degassing | bar |
| `f_O2` | Reference oxygen-fugacity buffer | `'IW'`, `'QFM'`, `'NiNiO'` |
| `buffer` | Offset from the reference buffer | log₁₀ units; e.g. `3` means buffer + 3 |
| `gas_behavior` | Gas treatment | `'ideal'` or `'real'` |
| `Degassing_style` | Degassing formulation | `'Intrusive'` or `'Extrusive'` |
| `Partitioning` | Calculate rock-to-melt volatile partitioning | `True` or `False` |
| `melting_style` | Partial-melting formulation | `'batch_melting'` or `'fractional_melting'` |
| `F_melt` | Fraction of source rock melted | Dimensionless; use `0 < F_melt <= 1` |

String settings are case-sensitive. In the solver functions, the buffer name is passed as `fugacity`; the variable called `f_O2` in `Set_up.py` contains this name rather than a numerical oxygen fugacity.

### Melt composition

Specify oxide abundances as **weight fractions**, using `mSiO2`, `mTiO2`, `mAl2O3`, `mFeO`, `mMgO`, `mCaO`, `mNa2O`, `mK2O`, and `mMnO`. For example, 48.3 wt% SiO₂ is entered as `mSiO2 = 0.483`.

The solver redistributes the supplied iron among Fe, FeO, and Fe₂O₃ and normalizes the resulting major-component composition internally. Basalt, alkali basalt, and rhyolite examples are included in `Set_up.py`.

### Volatile inventory

Volatile inputs are also **weight fractions**:

| Variable | Meaning when `Partitioning = False` | Default |
| --- | --- | --- |
| `mH2O` | Melt H₂O-equivalent content | `1000e-6` = 1000 ppm |
| `mCO2` | Melt CO₂-equivalent content | `500e-6` = 500 ppm |
| `mS` | Melt elemental sulfur content | `1000e-6` = 1000 ppm |
| `mN` | Melt elemental nitrogen content | `100e-6` = 100 ppm |

For conversion, `weight fraction = wt% / 100 = ppm / 1e6`.

When `Partitioning = True`, the volatile inputs are passed to the rock-to-melt partitioning routine. **The current implementation passes `mCO2` to the helper's elemental-carbon input, `m_C_mantle`.** In this mode, supply an elemental source-carbon weight fraction in that variable. The helper applies a redox-dependent carbon constraint and converts its result to melt CO₂. This convention differs from the CO₂-equivalent input used when partitioning is disabled.

## Degassing formulations

### Intrusive degassing

```python
Degassing_style = 'Intrusive'
```

The intrusive solver couples gas equilibrium, volatile mass balances, and solubilities determined by the vapor phase in equilibrium with the melt. Gas composition, dissolved volatile contents, and the gas-phase parameter `alphaG` are solved together.

The ten returned gas abundances are normalized to sum to one when gas species are present. They describe the **composition of the gas phase**. `alphaG` is the gas-phase mole-fraction parameter used alongside the melt fraction `1 - alphaG` in the solver's mass balances; it is not a gas species abundance or a gas mass fraction.

### Extrusive degassing

```python
Degassing_style = 'Extrusive'
```

The extrusive solver uses specified atmospheric partial pressures to determine solubilities. Set `PH2O`, `PH2`, `PCO2`, `PCO`, `PCH4`, `PH2S`, `PSO2`, `PS2`, `PNH3`, and `PN2` in **bar** in `Set_up.py`. The setup calculates total pressure from their sum.

The returned gas abundances are **not normalized to sum to one**. They are expressed on the solver's melt-based molar inventory scale, after its dissolution corrections. Normalize them separately only when a relative composition of the released gas is required. This solver does not return `alphaG` or update the atmosphere over time.

The current function signature accepts atmospheric partial pressures through `PN2`. `PNH3` contributes to total pressure in `Set_up.py` but is not passed separately to the solver.

**The zero-atmosphere settings in the supplied setup are not a working extrusive example.** See the implementation notes below before using this mode.

## Solubility options

Select the parameterizations with `Sol_H2O`, `Sol_CO2`, `Sol_S2`, `Sol_S6`, and `Sol_N`. These correspond to `H_law`, `C_law`, `S2_law`, `S6_law`, and `N_law` in the solver functions.

The following names follow the annotations in the source code:

| Setting | Supported keys | Parameterizations |
| --- | --- | --- |
| `Sol_H2O`, `Sol_CO2` | `'IM'` | Iacono Marziano (2012) |
| `Sol_H2O`, `Sol_CO2` | `'PW_basalt'`, `'PW_rhyolite'` | Parfitt and Wilson (2008) |
| `Sol_H2O`, `Sol_CO2` | `'LI'` | Lichtenberg (2021) |
| `Sol_H2O`, `Sol_CO2` | `'BO'` | Bower (2019) |
| `Sol_H2O`, `Sol_CO2` | `'DI'` | Dixon (1995, 1997) |
| `Sol_H2O`, `Sol_CO2` | `'BRA_basalt'`, `'BRA_rhyolite'` | Brachmann (2026) fits |
| `Sol_S2` | `'GAI'`, `'BOU'`, `'ONE'` | Gaillard (2026); Boulliung and Wood (2023); O'Neill (2021) |
| `Sol_S6` | `'SCA'`, `'BOU'`, `'ONE'` | Scaillet and Pichavant (2009); Boulliung and Wood (2022); O'Neill and Mavrogenes (2022) |
| `Sol_N` | `'BER'` | Bernadou (2021) |
| `Sol_N` | `'LIB_1'`, `'LIB_2'`, `'LIB_3'` | Libourel (2003); Libourel (2025) without/with NBO dependence |
| `Sol_N` | `'BRA'` | Brachmann (2026) fit |
| `Sol_N` | `'DAS'` | Additional nitrogen parameterization implemented in the solvers |

The defaults are `'BRA_basalt'` for H₂O and CO₂, `'GAI'` for sulfide, `'BOU'` for sulfate, and `'BRA'` for nitrogen. Sulfide or sulfate solubility is selected internally according to the redox buffer and offset.

Choose laws appropriate to the melt composition and their calibration conditions. A selectable key does not establish validity across all pressures, temperatures, or compositions.

## Using the solver from Python

Run the following example from the repository directory, or save it there as a separate script:

```python
from solve_gases_CHONS import solve_gases_CHONS

result = solve_gases_CHONS(
    P_melt=100, T_melt=1675, fugacity='IW', buffer=3,
    mH2O_init=1000e-6, mCO2_init=500e-6, mS_init=1000e-6, mN_init=100e-6,
    H_law='BRA_basalt', C_law='BRA_basalt', S2_law='GAI', S6_law='BOU', N_law='BRA',
    mSiO2=0.483, mTiO2=0.013, mAl2O3=0.179, mFeO=0.088, mMgO=0.081,
    mCaO=0.103, mNa2O=0.024, mK2O=0.002, mMnO=0,
    Partitioning=False, melting_style='batch_melting', F_melt=0.25, gas_behavior='ideal'
)
names = (
    'H2', 'H2O', 'CH4', 'CO', 'CO2', 'S2', 'SO2', 'H2S', 'N2', 'NH3',
    'H2O_melt', 'CO2_melt', 'S_melt', 'N_melt', 'alphaG'
)
results = dict(zip(names, result))
for name, value in results.items():
    print(f'{name}: {value:.6e}')
```

The intrusive solver returns 15 values in the order shown above. The extrusive function, `solve_gases_CHONS_ext` in `solve_gases_CHONS_extrusive.py`, returns the same first 14 values and omits `alphaG`.

`XH2O_melt`, `XCO2_melt`, `XS_melt`, and `XN_melt` are **melt mole-fraction outputs**, not weight fractions or ppm. Sulfur and nitrogen use elemental S and N accounting. Converting these outputs to mass units requires the molecular weights and the melt normalization used by the solver.

## Repository files

| File | Purpose |
| --- | --- |
| [`Set_up.py`](Set_up.py) | User settings and example solver calls |
| [`solve_gases_CHONS.py`](solve_gases_CHONS.py) | Coupled intrusive degassing solver |
| [`solve_gases_CHONS_extrusive.py`](solve_gases_CHONS_extrusive.py) | Extrusive degassing solver using atmospheric partial pressures |
| [`initial_guess.py`](initial_guess.py) | Initial estimates for the intrusive solver |
| [`initial_guess_ext.py`](initial_guess_ext.py) | Initial estimates for the extrusive solver |
| [`func_CHNS_part.py`](func_CHNS_part.py) | Volatile partitioning during partial melting |
| [`func_Fe_FeO_Fe2O3.py`](func_Fe_FeO_Fe2O3.py) | Iron redox distribution |
| [`Gibbs_E.py`](Gibbs_E.py) | Thermodynamic functions used to calculate equilibrium constants |
| [`FugacityC_ALL.py`](FugacityC_ALL.py) | Approximate fugacity coefficients |
| [`Molecular_weights.py`](Molecular_weights.py) | Molecular and atomic weights |

## Implementation notes and troubleshooting

- **Initial-guess pressure range:** Both initial-guess routines interpolate over **10⁻⁶–10⁴ bar**. Calls outside this interval raise an interpolation error. This is a numerical input range, not a scientific validity range.
- **Zero atmospheric pressure:** `Set_up.py` replaces a zero total pressure with `1e-10` bar, which is below that interpolation range. Several extrusive solubility expressions also evaluate logarithms of atmospheric partial pressures directly. Zero inputs can therefore generate infinities or NaNs; the zero-atmosphere case requires code-level handling.
- **Approximate real-gas treatment:** `FugacityC_ALL.py` describes its coefficients as approximations and warns against pressures above **5000 bar**. N₂ and NH₃ do not receive explicit fugacity-coefficient corrections in the main solvers. The Dixon solubility options call fugacity calculations even when `gas_behavior = 'ideal'`.
- **Convergence:** The solvers use SciPy nonlinear least squares and raise `RuntimeError` when the optimizer reports failure. Optimizer success alone does not verify the physical solution; check element balances, residuals, and the applicability of the chosen laws when evaluating new conditions. The intrusive gas outputs are normalized after solving, so their sum alone is not a convergence test.

## License

MORDOR is distributed under the **MIT License**. See [`LICENSE`](LICENSE) for the full terms.
