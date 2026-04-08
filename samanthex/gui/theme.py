"""GUI theme and styling for SAMANTHEX.

Defines a professional dark theme using ttk styles, optimized for
laboratory environments with reduced eye strain. All colors, fonts,
and widget styles are centralized here for consistency.
"""

import tkinter as tk
from tkinter import ttk


# --- Color Palette ---
class Colors:
    """Dark theme color palette."""
    BG_PRIMARY = "#1e1e2e"       # Main background
    BG_SECONDARY = "#2b2b3d"     # Panel/frame background
    BG_TERTIARY = "#363649"      # Input fields, cards
    BG_HOVER = "#414158"         # Hover state

    FG_PRIMARY = "#cdd6f4"       # Main text
    FG_SECONDARY = "#a6adc8"     # Secondary/muted text
    FG_DISABLED = "#585b70"      # Disabled text

    ACCENT = "#89b4fa"           # Primary accent (blue)
    ACCENT_HOVER = "#74c7ec"     # Accent hover
    SUCCESS = "#a6e3a1"          # Success / good quality
    WARNING = "#f9e2af"          # Warning / medium quality
    ERROR = "#f38ba8"            # Error / low quality
    INFO = "#89dceb"             # Informational

    BORDER = "#45475a"           # Borders
    SELECTION = "#313244"        # Selected items

    # Quality colors for visualization
    Q_HIGH = "#a6e3a1"           # Q >= 30 (green)
    Q_MEDIUM = "#f9e2af"         # Q >= 20 (yellow)
    Q_LOW = "#f38ba8"            # Q < 20 (red)

    # Button colors
    BTN_PRIMARY_BG = "#89b4fa"
    BTN_PRIMARY_FG = "#1e1e2e"
    BTN_SUCCESS_BG = "#a6e3a1"
    BTN_SUCCESS_FG = "#1e1e2e"
    BTN_DANGER_BG = "#f38ba8"
    BTN_DANGER_FG = "#1e1e2e"


# --- Fonts ---
class Fonts:
    """Font definitions."""
    FAMILY = "Segoe UI"
    MONO = "Consolas"

    TITLE = (FAMILY, 16, "bold")
    SUBTITLE = (FAMILY, 12, "bold")
    HEADING = (FAMILY, 11, "bold")
    BODY = (FAMILY, 10)
    SMALL = (FAMILY, 9)
    BUTTON = (FAMILY, 10, "bold")

    CODE = (MONO, 10)
    CODE_SMALL = (MONO, 9)


def configure_styles(root: tk.Tk) -> ttk.Style:
    """Configure ttk styles for the dark theme.

    Args:
        root: The root Tk window.

    Returns:
        Configured ttk.Style instance.
    """
    style = ttk.Style(root)

    # Use 'clam' as base — most customizable cross-platform theme
    style.theme_use("clam")

    # --- General ---
    style.configure(".", background=Colors.BG_PRIMARY, foreground=Colors.FG_PRIMARY,
                    font=Fonts.BODY, borderwidth=0, relief="flat")

    # --- Frame ---
    style.configure("TFrame", background=Colors.BG_PRIMARY)
    style.configure("Card.TFrame", background=Colors.BG_SECONDARY)

    # --- Label ---
    style.configure("TLabel", background=Colors.BG_PRIMARY,
                    foreground=Colors.FG_PRIMARY, font=Fonts.BODY)
    style.configure("Title.TLabel", font=Fonts.TITLE, foreground=Colors.ACCENT)
    style.configure("Subtitle.TLabel", font=Fonts.SUBTITLE)
    style.configure("Heading.TLabel", font=Fonts.HEADING)
    style.configure("Small.TLabel", font=Fonts.SMALL, foreground=Colors.FG_SECONDARY)
    style.configure("Success.TLabel", foreground=Colors.SUCCESS)
    style.configure("Warning.TLabel", foreground=Colors.WARNING)
    style.configure("Error.TLabel", foreground=Colors.ERROR)
    style.configure("Status.TLabel", background=Colors.BG_SECONDARY,
                    foreground=Colors.FG_SECONDARY, font=Fonts.SMALL, padding=(8, 4))

    # --- LabelFrame ---
    style.configure("TLabelframe", background=Colors.BG_PRIMARY,
                    foreground=Colors.FG_PRIMARY, bordercolor=Colors.BORDER,
                    relief="solid", borderwidth=1)
    style.configure("TLabelframe.Label", background=Colors.BG_PRIMARY,
                    foreground=Colors.ACCENT, font=Fonts.HEADING)

    # --- Button ---
    style.configure("TButton", background=Colors.BG_TERTIARY,
                    foreground=Colors.FG_PRIMARY, font=Fonts.BUTTON,
                    padding=(12, 6), borderwidth=0)
    style.map("TButton",
              background=[("active", Colors.BG_HOVER), ("disabled", Colors.BG_SECONDARY)],
              foreground=[("disabled", Colors.FG_DISABLED)])

    style.configure("Primary.TButton", background=Colors.BTN_PRIMARY_BG,
                    foreground=Colors.BTN_PRIMARY_FG)
    style.map("Primary.TButton",
              background=[("active", Colors.ACCENT_HOVER), ("disabled", Colors.BG_SECONDARY)],
              foreground=[("disabled", Colors.FG_DISABLED)])

    style.configure("Success.TButton", background=Colors.BTN_SUCCESS_BG,
                    foreground=Colors.BTN_SUCCESS_FG)
    style.map("Success.TButton",
              background=[("active", "#94d990"), ("disabled", Colors.BG_SECONDARY)],
              foreground=[("disabled", Colors.FG_DISABLED)])

    # --- Entry ---
    style.configure("TEntry", fieldbackground=Colors.BG_TERTIARY,
                    foreground=Colors.FG_PRIMARY, insertcolor=Colors.FG_PRIMARY,
                    bordercolor=Colors.BORDER, lightcolor=Colors.BG_TERTIARY,
                    darkcolor=Colors.BG_TERTIARY, borderwidth=1, padding=(6, 4))
    style.map("TEntry",
              fieldbackground=[("focus", Colors.BG_HOVER)],
              bordercolor=[("focus", Colors.ACCENT)])

    # --- Combobox ---
    style.configure("TCombobox", fieldbackground=Colors.BG_TERTIARY,
                    foreground=Colors.FG_PRIMARY, background=Colors.BG_TERTIARY,
                    arrowcolor=Colors.FG_PRIMARY, bordercolor=Colors.BORDER,
                    lightcolor=Colors.BG_TERTIARY, darkcolor=Colors.BG_TERTIARY,
                    borderwidth=1, padding=(6, 4))
    style.map("TCombobox",
              fieldbackground=[("readonly", Colors.BG_TERTIARY)],
              foreground=[("readonly", Colors.FG_PRIMARY)],
              bordercolor=[("focus", Colors.ACCENT)])

    # --- Notebook ---
    style.configure("TNotebook", background=Colors.BG_PRIMARY,
                    borderwidth=0)
    style.configure("TNotebook.Tab", background=Colors.BG_SECONDARY,
                    foreground=Colors.FG_SECONDARY, font=Fonts.BODY,
                    padding=(16, 6), borderwidth=0)
    style.map("TNotebook.Tab",
              background=[("selected", Colors.BG_PRIMARY)],
              foreground=[("selected", Colors.ACCENT)])

    # --- Progressbar ---
    style.configure("TProgressbar", background=Colors.ACCENT,
                    troughcolor=Colors.BG_SECONDARY, borderwidth=0,
                    thickness=6)
    style.configure("Green.Horizontal.TProgressbar", background=Colors.SUCCESS,
                    troughcolor=Colors.BG_SECONDARY, borderwidth=0)

    # --- Scrollbar ---
    style.configure("TScrollbar", background=Colors.BG_SECONDARY,
                    troughcolor=Colors.BG_PRIMARY, borderwidth=0,
                    arrowcolor=Colors.FG_SECONDARY)
    style.map("TScrollbar", background=[("active", Colors.BG_HOVER)])

    # --- Checkbutton ---
    style.configure("TCheckbutton", background=Colors.BG_PRIMARY,
                    foreground=Colors.FG_PRIMARY, font=Fonts.BODY)
    style.map("TCheckbutton",
              background=[("active", Colors.BG_PRIMARY)])

    # --- Separator ---
    style.configure("TSeparator", background=Colors.BORDER)

    return style
