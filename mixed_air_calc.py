"""
Author: Alex Kalmbach

Description: This program determines the state of mixed air, given the properties of two separate air streams.
Each airstream is represented as an object of the AirStream class. So OOP is being implemented in this code!
Adiabatic mixing is assumed (i.e. no heat transfer from surroundings). The mixed air state is calculated from
conservation laws -- specifically, the conservation of mass (for air and water vapor) and conservation of energy (enthalpy).
The user can choose to select the airflow values in actual CFM (ACFM) or standard CFM (SCFM).
All psychrometric mass and energy calculations are based on the ACFM value.

The CoolProp library is implemented to obtain humid air properties.

References:
    [1] Colmac Coil Useful Equations When Sizing Coils by Carl E. Berg
    [2] Price Industries Engineer's HVAC Handbook 2023 (Chapter 5 - Psychrometrics)

Date: Version 1.0 Completed 02/03/2025
"""

import tkinter as tk
from tkinter import ttk
import CoolProp.CoolProp as CP
from pint import UnitRegistry
import tkinter.messagebox as mbox

ureg = UnitRegistry()
Q_ = ureg.Quantity

STD_AIR_DENSITY = 0.075  # lb/ft^3

def scfm_to_acfm(scfm: float, rho_act: float) -> float:
    """
    Convert SCFM to ACFM using density. From Colmac article [1]
      ACFM = SCFM * (rho_std / rho_act)
    """
    return scfm * (STD_AIR_DENSITY / rho_act)

def acfm_to_scfm(acfm: float, rho_act: float) -> float:
    """
    Convert ACFM to SCFM (inverse). From Colmac article [1]
      SCFM = ACFM * (rho_act / rho_std)
    """
    return acfm * (rho_act / STD_AIR_DENSITY)

def convert_to_si(property_name: str, value: float) -> float:
    conversions = {
        "Tdb": Q_(value, "degF").to("kelvin").magnitude,
        "Twb": Q_(value, "degF").to("kelvin").magnitude,
        "Tdp": Q_(value, "degF").to("kelvin").magnitude,
        "P":   Q_(value, "psi").to("pascal").magnitude,
        "RH":  value / 100.0,  # Convert % to fraction
        "W":   Q_(value, "gr/lb").to("kg/kg").magnitude,
        "H":   Q_(value, "Btu/lb").to("J/kg").magnitude,
    }
    return conversions.get(property_name, value)

def convert_from_si(property_name: str, value: float) -> float:
    conversions = {
        "Tdb":    Q_(value, "kelvin").to("degF").magnitude,
        "Twb":    Q_(value, "kelvin").to("degF").magnitude,
        "Tdp":    Q_(value, "kelvin").to("degF").magnitude,
        "P":      Q_(value, "pascal").to("psi").magnitude,
        "RH":     value * 100.0,  # fraction -> %
        "W":      Q_(value, "kg/kg").to("gr/lb").magnitude,
        "V":      Q_(value, "m^3/kg").to("ft^3/lb").magnitude,
        "H":      Q_(value, "J/kg").to("Btu/lb").magnitude,
        "D":      Q_(value, "kg/m^3").to("lb/ft^3").magnitude,
    }
    return conversions.get(property_name, value)

class AirStream:
    def __init__(
        self,
        prop1_name: str,
        prop1_value: float,
        prop2_name: str,
        prop2_value: float,
        pressure_psi: float,
        flow_type: str,   # "ACFM" or "SCFM"
        flow_value: float
    ):
        self.pressure_psi = pressure_psi
        self.pressure_si = convert_to_si("P", pressure_psi)

        # Convert the input psychrometric properties to SI
        self.prop1_name = prop1_name
        self.prop1_value_si = convert_to_si(prop1_name, prop1_value)
        self.prop2_name = prop2_name
        self.prop2_value_si = convert_to_si(prop2_name, prop2_value)

        # Calculate derived properties in SI, then convert to IP
        self.properties = self.calculate_properties()

        # Get the actual humid-air density [lb/ft^3]
        rho_act = self.properties["D"]

        # Figure out ACFM/SCFM based on user selection
        flow_type = flow_type.strip().upper()
        if flow_type == "SCFM":
            self.scfm = flow_value
            self.acfm = scfm_to_acfm(self.scfm, rho_act)
        elif flow_type == "ACFM":
            self.acfm = flow_value
            self.scfm = acfm_to_scfm(self.acfm, rho_act)
        else:
            raise ValueError("flow_type must be either 'ACFM' or 'SCFM'.")

        # For mass/energy balances, we use ACFM flow
        # air_density is lb of DRY air per ft^3
        air_density = 1.0 / self.properties["V"]  # lb/ft^3 of dry air
        # lb/min of DRY air
        self.air_mass_flow_rate = air_density * self.acfm
        # lb/min of water vapor
        self.water_mass_flow_rate = (
            self.air_mass_flow_rate * (self.properties["W"] / 7000.0)
        )

    def calculate_properties(self) -> dict:
        keys = ["Tdb", "Twb", "Tdp", "RH", "W", "V", "H"]
        prop_dict_si = {}
        try:
            for k in keys:
                prop_dict_si[k] = CP.HAPropsSI(
                    k,
                    self.prop1_name, self.prop1_value_si,
                    self.prop2_name, self.prop2_value_si,
                    "P", self.pressure_si
                )
        except Exception as e:
            # Trim the error if "inputs were" is present
            e_str = str(e)
            idx = e_str.lower().find("inputs were")
            if idx != -1:
                e_str = e_str[:idx].strip(" :\n\r\t")
            raise ValueError(f"CoolProp error: {e_str}")

        # Add density in SI, convert to IP
        prop_dict_si["D"] = 1.0 / prop_dict_si["V"]
        prop_dict_ip = {k: convert_from_si(k, v) for k, v in prop_dict_si.items()}
        return prop_dict_ip

    @staticmethod
    def calculate_pressure_from_elevation(elev_ft: float) -> float:
        """
        Eqn. 5.17 in the Price HVAC Handbook [2]. Converts elevation to pressure in psia.
        """
        return 14.696 * (1 - 6.8754e-6 * elev_ft) ** 5.2559

class MixedAirApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Adiabatic Mixing of Airstreams")

        # 1) Define styles for a more professional appearance
        style = ttk.Style(self)
        style.configure("TLabel",     font=("Helvetica", 10))
        style.configure("TEntry",     font=("Helvetica", 10))
        style.configure("TButton",    font=("Helvetica", 10))
        style.configure("TCombobox",  font=("Helvetica", 10))
        style.configure("TLabelframe", padding="10 10 10 10")
        style.configure("TLabelframe.Label", font=("Helvetica", 10, "bold"))
        style.configure("Title.TLabel", font=("Helvetica", 16, "bold"))

        # Use grid for the main window
        self.columnconfigure(0, weight=1)

        # ----------------------------
        #  Top Row: Title & Info Button
        # ----------------------------
        top_row_frame = ttk.Frame(self, padding="10 10 10 0")
        top_row_frame.grid(row=0, column=0, sticky="ew")
        top_row_frame.columnconfigure(0, weight=1)

        self.title_label = ttk.Label(
            top_row_frame,
            text="Adiabatic Mixing of Airstreams",
            style="Title.TLabel"
        )
        self.title_label.grid(row=0, column=0, sticky="w")

        # Info button on the right
        info_button = ttk.Button(top_row_frame, text="Info.", command=self._show_info)
        info_button.grid(row=0, column=1, sticky="e")

        # ----------------------------
        #  Row 1: Error Label
        # ----------------------------
        self.error_label = ttk.Label(self, text="", foreground="red", padding="0 5 0 0")
        self.error_label.grid(row=1, column=0, sticky="ew")

        # ----------------------------
        #  Row 2: Pressure/Elevation Controls
        # ----------------------------
        self.top_frame = ttk.Frame(self, padding="10 10 10 10")
        self.top_frame.grid(row=2, column=0, sticky="ew")

        self.pressure_mode_var = tk.StringVar(value="psi")
        pressure_mode_cb = ttk.Combobox(
            self.top_frame,
            textvariable=self.pressure_mode_var,
            values=["psi", "ft"],  # psia or elevation in ft
            width=5,
            state="readonly"
        )
        pressure_mode_cb.grid(row=0, column=0, padx=(0, 5), pady=5)
        pressure_mode_cb.bind("<<ComboboxSelected>>", self._on_input_change)

        self.pressure_label_var = tk.StringVar(value="Pressure (psia):")
        self.pressure_label = ttk.Label(self.top_frame, textvariable=self.pressure_label_var)
        self.pressure_label.grid(row=0, column=1, sticky="e")

        self.pressure_value_var = tk.StringVar(value="14.7")
        def only_numeric_input(P):
            if P in ["", "-", ".", "-."]:
                return True
            try:
                float(P)
                return True
            except ValueError:
                return False

        self.vcmd = (self.register(only_numeric_input), "%P")

        pressure_entry = ttk.Entry(
            self.top_frame,
            textvariable=self.pressure_value_var,
            width=10,
            validate="key",
            validatecommand=self.vcmd
        )
        pressure_entry.grid(row=0, column=2, padx=(5, 10))
        self.pressure_value_var.trace_add("write", self._on_input_change)

        # ----------------------------
        #  Row 3: Airstreams + Mixed
        # ----------------------------
        self.main_content_frame = ttk.Frame(self, padding="5 10 5 10")
        self.main_content_frame.grid(row=3, column=0, sticky="nsew")
        self.main_content_frame.columnconfigure(0, weight=1)
        self.main_content_frame.columnconfigure(1, weight=1)
        self.main_content_frame.columnconfigure(2, weight=1)

        # AIRSTREAM #1
        self.airstream1_frame = ttk.LabelFrame(self.main_content_frame, text="AirStream #1")
        self.airstream1_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        row_idx = 0
        ttk.Label(self.airstream1_frame, text="Tdb (°F):").grid(row=row_idx, column=0, sticky="e")
        self.as1_tdb_var = tk.StringVar(value="75.0")
        as1_tdb_entry = ttk.Entry(
            self.airstream1_frame,
            textvariable=self.as1_tdb_var,
            width=8,
            validate="key",
            validatecommand=self.vcmd
        )
        as1_tdb_entry.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as1_tdb_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream1_frame, text="2nd Property:").grid(row=row_idx, column=0, sticky="e")
        self.prop_choices = [
            ("RH (%)",    "RH"),
            ("Twb (°F)",  "Twb"),
            ("Tdp (°F)",  "Tdp"),
            ("W (gr/lb)", "W")
        ]
        self.prop_map = { disp: key for (disp, key) in self.prop_choices }

        self.as1_prop2_display_var = tk.StringVar(value="RH (%)")
        as1_prop2_name_cb = ttk.Combobox(
            self.airstream1_frame,
            textvariable=self.as1_prop2_display_var,
            values=[p[0] for p in self.prop_choices],
            width=10,
            state="readonly"
        )
        as1_prop2_name_cb.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as1_prop2_display_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream1_frame, text="Value:").grid(row=row_idx, column=0, sticky="e")
        self.as1_prop2_value_var = tk.StringVar(value="55.0")
        as1_prop2_value_entry = ttk.Entry(
            self.airstream1_frame,
            textvariable=self.as1_prop2_value_var,
            width=8,
            validate="key",
            validatecommand=self.vcmd
        )
        as1_prop2_value_entry.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as1_prop2_value_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream1_frame, text="Flow Type:").grid(row=row_idx, column=0, sticky="e")
        self.as1_flow_type_var = tk.StringVar(value="SCFM")
        as1_flow_type_cb = ttk.Combobox(
            self.airstream1_frame,
            textvariable=self.as1_flow_type_var,
            values=["SCFM", "ACFM"],
            width=6,
            state="readonly"
        )
        as1_flow_type_cb.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as1_flow_type_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream1_frame, text="Flow (cfm):").grid(row=row_idx, column=0, sticky="e")
        self.as1_flow_value_var = tk.StringVar(value="1500")
        as1_flow_value_entry = ttk.Entry(
            self.airstream1_frame,
            textvariable=self.as1_flow_value_var,
            width=8,
            validate="key",
            validatecommand=self.vcmd
        )
        as1_flow_value_entry.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as1_flow_value_var.trace_add("write", self._on_input_change)

        row_idx += 1
        self.as1_output_frame = ttk.Frame(self.airstream1_frame)
        self.as1_output_frame.grid(row=row_idx, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        # Updated props_list to include mass flows
        props_list = ["Tdb", "Twb", "Tdp", "RH", "W", "H", "D", "ACFM", "SCFM", "mAir", "mWater"]
        self.as1_labels = {}
        r, c = 0, 0
        for prop in props_list:
            self.as1_labels[prop] = ttk.Label(self.as1_output_frame, text=f"{prop}: ---")
            self.as1_labels[prop].grid(row=r, column=c, padx=4, pady=2, sticky="w")
            c += 1
            if c == 3:
                c = 0
                r += 1

        # AIRSTREAM #2
        self.airstream2_frame = ttk.LabelFrame(self.main_content_frame, text="AirStream #2")
        self.airstream2_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

        row_idx = 0
        ttk.Label(self.airstream2_frame, text="Tdb (°F):").grid(row=row_idx, column=0, sticky="e")
        self.as2_tdb_var = tk.StringVar(value="95.0")
        as2_tdb_entry = ttk.Entry(
            self.airstream2_frame,
            textvariable=self.as2_tdb_var,
            width=8,
            validate="key",
            validatecommand=self.vcmd
        )
        as2_tdb_entry.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as2_tdb_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream2_frame, text="2nd Property:").grid(row=row_idx, column=0, sticky="e")
        self.as2_prop2_display_var = tk.StringVar(value="Twb (°F)")
        as2_prop2_name_cb = ttk.Combobox(
            self.airstream2_frame,
            textvariable=self.as2_prop2_display_var,
            values=[p[0] for p in self.prop_choices],
            width=10,
            state="readonly"
        )
        as2_prop2_name_cb.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as2_prop2_display_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream2_frame, text="Value:").grid(row=row_idx, column=0, sticky="e")
        self.as2_prop2_value_var = tk.StringVar(value="75.0")
        as2_prop2_value_entry = ttk.Entry(
            self.airstream2_frame,
            textvariable=self.as2_prop2_value_var,
            width=8,
            validate="key",
            validatecommand=self.vcmd
        )
        as2_prop2_value_entry.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as2_prop2_value_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream2_frame, text="Flow Type:").grid(row=row_idx, column=0, sticky="e")
        self.as2_flow_type_var = tk.StringVar(value="SCFM")
        as2_flow_type_cb = ttk.Combobox(
            self.airstream2_frame,
            textvariable=self.as2_flow_type_var,
            values=["SCFM", "ACFM"],
            width=6,
            state="readonly"
        )
        as2_flow_type_cb.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as2_flow_type_var.trace_add("write", self._on_input_change)

        row_idx += 1
        ttk.Label(self.airstream2_frame, text="Flow (cfm):").grid(row=row_idx, column=0, sticky="e")
        self.as2_flow_value_var = tk.StringVar(value="500")
        as2_flow_value_entry = ttk.Entry(
            self.airstream2_frame,
            textvariable=self.as2_flow_value_var,
            width=8,
            validate="key",
            validatecommand=self.vcmd
        )
        as2_flow_value_entry.grid(row=row_idx, column=1, padx=5, pady=2)
        self.as2_flow_value_var.trace_add("write", self._on_input_change)

        row_idx += 1
        self.as2_output_frame = ttk.Frame(self.airstream2_frame)
        self.as2_output_frame.grid(row=row_idx, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.as2_labels = {}
        r, c = 0, 0
        for prop in props_list:
            self.as2_labels[prop] = ttk.Label(self.as2_output_frame, text=f"{prop}: ---")
            self.as2_labels[prop].grid(row=r, column=c, padx=4, pady=2, sticky="w")
            c += 1
            if c == 3:
                c = 0
                r += 1

        # MIXED STREAM
        self.mixed_frame = ttk.LabelFrame(self.main_content_frame, text="Mixed Air Stream")
        self.mixed_frame.grid(row=0, column=2, padx=5, pady=5, sticky="nsew")

        self.mixed_labels_frame = ttk.Frame(self.mixed_frame)
        self.mixed_labels_frame.pack(pady=5, fill="x")

        self.mixed_labels = {}
        for prop in props_list:
            self.mixed_labels[prop] = ttk.Label(self.mixed_labels_frame, text=f"{prop}: ---")
            self.mixed_labels[prop].pack(anchor="w", padx=4, pady=2)

        # ----------------------------
        #  Row 4: Author Credit
        # ----------------------------
        self.author_label = ttk.Label(
            self,
            text="Developed in Python by Alex Kalmbach",
            foreground="gray"
        )
        self.author_label.grid(row=4, column=0, sticky="e", padx=5, pady=5)

        # Initial calculation
        self._on_input_change()

    def _show_info(self):
        mbox.showinfo(
            "Information",
            "This program calculates the adiabatic mixing of two airstreams "
            "by applying the laws of conservation of mass and energy.\n\n"
            f"The user can enter airflow at standard cfm (SCFM), "
            f"which has a density of {STD_AIR_DENSITY} lb/ft³, "
            "or actual cfm (ACFM), which has a density based on the air at user-specified conditions.\n\n"
            "The CoolProp library is used to calculate psychrometric properties of air.\n\n"
            "Note: Enthalpy values have a zero Btu/lb reference based on SI units, which may not align with IP based psych charts."
        )

    def _on_input_change(self, *args):
        self.error_label.config(text="")
        pressure_str = self.pressure_value_var.get().strip()
        if not pressure_str:
            return
        try:
            if self.pressure_mode_var.get() == "psi":
                p_psi = float(pressure_str)
            else:
                elevation_ft = float(pressure_str)
                p_psi = AirStream.calculate_pressure_from_elevation(elevation_ft)
        except ValueError as e:
            self.error_label.config(text=f"Invalid pressure/elevation: {e}")
            return

        label_txt = "Pressure (psia):" if (self.pressure_mode_var.get() == "psi") else "Elevation (ft):"
        self.pressure_label_var.set(label_txt)

        # Airstream #1
        try:
            as1_tdb = float(self.as1_tdb_var.get())
            as1_prop2_name = self.prop_map[self.as1_prop2_display_var.get()]
            as1_prop2_val = float(self.as1_prop2_value_var.get())
            as1_flow_type = self.as1_flow_type_var.get()
            as1_flow_val = float(self.as1_flow_value_var.get())
            if as1_prop2_name.lower() == "rh" and not (0 <= as1_prop2_val <= 100):
                raise ValueError("RH must be between 0 and 100%")

            a1 = AirStream(
                prop1_name="Tdb",
                prop1_value=as1_tdb,
                prop2_name=as1_prop2_name,
                prop2_value=as1_prop2_val,
                pressure_psi=p_psi,
                flow_type=as1_flow_type,
                flow_value=as1_flow_val
            )
            self._update_airstream_labels(a1, self.as1_labels)
        except ValueError as e:
            self.error_label.config(text=f"Airstream #1 invalid: {e}")
            return

        # Airstream #2
        try:
            as2_tdb = float(self.as2_tdb_var.get())
            as2_prop2_name = self.prop_map[self.as2_prop2_display_var.get()]
            as2_prop2_val = float(self.as2_prop2_value_var.get())
            as2_flow_type = self.as2_flow_type_var.get()
            as2_flow_val = float(self.as2_flow_value_var.get())

            if as2_prop2_name.lower() == "rh" and not (0 <= as2_prop2_val <= 100):
                raise ValueError("RH must be between 0 and 100%")

            a2 = AirStream(
                prop1_name="Tdb",
                prop1_value=as2_tdb,
                prop2_name=as2_prop2_name,
                prop2_value=as2_prop2_val,
                pressure_psi=p_psi,
                flow_type=as2_flow_type,
                flow_value=as2_flow_val
            )
            self._update_airstream_labels(a2, self.as2_labels)
        except ValueError as e:
            self.error_label.config(text=f"Airstream #2 invalid: {e}")
            return

        # Mixed Stream
        try:
            mixed = self._calculate_mixed_stream(a1, a2, p_psi)
            self._update_mixed_labels(mixed)
        except ValueError as e:
            self.error_label.config(text=f"Mixed stream error: {e}")
            return

    def _update_airstream_labels(self, airstream_obj, label_dict):
        props = airstream_obj.properties
        label_dict["Tdb"].config(text=f"Tdb: {props['Tdb']:.1f} °F")
        label_dict["Twb"].config(text=f"Twb: {props['Twb']:.1f} °F")
        label_dict["Tdp"].config(text=f"Tdp: {props['Tdp']:.1f} °F")
        label_dict["RH"].config(text=f"RH: {props['RH']:.1f} %")
        label_dict["W"].config(text=f"W: {props['W']:.2f} gr/lb")
        label_dict["H"].config(text=f"H: {props['H']:.2f} Btu/lb")
        label_dict["D"].config(text=f"D: {props['D']:.4f} lb/ft³")
        label_dict["ACFM"].config(text=f"ACFM: {airstream_obj.acfm:.0f} cfm")
        label_dict["SCFM"].config(text=f"SCFM: {airstream_obj.scfm:.0f} cfm")

        # Display air mass flow and water mass flow
        label_dict["mAir"].config(
            text=f"Dry Air Mass Flow: {airstream_obj.air_mass_flow_rate:.2f} lb/min"
        )
        label_dict["mWater"].config(
            text=f"Water Vapor Mass Flow: {airstream_obj.water_mass_flow_rate:.2f} lb/min"
        )

    def _calculate_mixed_stream(self, a1, a2, p_psi):
        """
        This function is the heart of the program! The mass and energy balance are what this program is based on.
        By solving for enthalpy (H) and humidity ratio (W), the other properties can be obtained. From [2]:
        m1 + m2 = m3
        m1w1 + m2w2 = m3w3
        m1h1 + m2h2 = m3mh3
        """
        air_mass_flow_mixed = a1.air_mass_flow_rate + a2.air_mass_flow_rate
        water_flow_mixed = a1.water_mass_flow_rate + a2.water_mass_flow_rate
        W_mixed = (water_flow_mixed / air_mass_flow_mixed) * 7000.0

        H_mixed = (
            a1.air_mass_flow_rate * a1.properties["H"] +
            a2.air_mass_flow_rate * a2.properties["H"]
        ) / air_mass_flow_mixed

        temp_mixed = AirStream(
            prop1_name="W", prop1_value=W_mixed,
            prop2_name="H", prop2_value=H_mixed,
            pressure_psi=p_psi,
            flow_type="ACFM",
            flow_value=0.0
        )
        mixed_density = temp_mixed.properties["D"]
        mixed_acfm = air_mass_flow_mixed / mixed_density

        final_mixed = AirStream(
            prop1_name="W", prop1_value=W_mixed,
            prop2_name="H", prop2_value=H_mixed,
            pressure_psi=p_psi,
            flow_type="ACFM",
            flow_value=mixed_acfm
        )
        return final_mixed

    def _update_mixed_labels(self, mixed):
        props = mixed.properties
        self.mixed_labels["Tdb"].config(text=f"Tdb: {props['Tdb']:.1f} °F")
        self.mixed_labels["Twb"].config(text=f"Twb: {props['Twb']:.1f} °F")
        self.mixed_labels["Tdp"].config(text=f"Tdp: {props['Tdp']:.1f} °F")
        self.mixed_labels["RH"].config(text=f"RH: {props['RH']:.1f} %")
        self.mixed_labels["W"].config(text=f"W: {props['W']:.2f} gr/lb")
        self.mixed_labels["H"].config(text=f"H: {props['H']:.2f} Btu/lb")
        self.mixed_labels["D"].config(text=f"D: {props['D']:.4f} lb/ft³")
        self.mixed_labels["ACFM"].config(text=f"ACFM: {mixed.acfm:.0f} cfm")
        self.mixed_labels["SCFM"].config(text=f"SCFM: {mixed.scfm:.0f} cfm")

        # Display air mass flow and water mass flow for mixed stream
        self.mixed_labels["mAir"].config(
            text=f"Dry Air Mass Flow: {mixed.air_mass_flow_rate:.2f} lb/min"
        )
        self.mixed_labels["mWater"].config(
            text=f"Water Vapor Mass Flow: {mixed.water_mass_flow_rate:.2f} lb/min"
        )

if __name__ == "__main__":
    app = MixedAirApp()
    app.mainloop()
