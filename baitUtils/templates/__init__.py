"""Static assets (CSS and JavaScript) for the HTML reports."""

from importlib import resources


def load_template(name: str) -> str:
    """Return the text of a template file shipped with the package."""
    return resources.files(__package__).joinpath(name).read_text(encoding="utf-8")
