import os

try:
    from importlib.resources import files as resource_files
except ImportError:
    try:
        from importlib_resources import files as resource_files
    except ImportError:
        resource_files = None

INFERNAL_FOLDER = os.path.realpath(
    os.path.join(
        __file__,
        os.pardir,
        os.pardir,
        os.pardir,
        os.pardir,
        "vendor",
        "infernal",
    )
)