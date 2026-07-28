#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description="Compile the Coherent Thomson Scattering project and run the solver."
    )
    parser.add_argument(
        "config_file",
        nargs="?",
        default=None,
        help="Path to the optional configuration file."
    )
    args = parser.parse_args()

    # Determine paths
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    build_dir = project_root / "build"
    binary_path = project_root / "bin" / "coherent_thomson_solver"

    print(f"Project root: {project_root}")
    print(f"Build directory: {build_dir}")

    # Step 1: Configure using CMake
    print("\n--- Configuring project with CMake ---")
    configure_cmd = ["cmake", "-S", str(project_root), "-B", str(build_dir)]
    try:
        subprocess.run(configure_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: CMake configuration failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)
    except FileNotFoundError:
        print("Error: cmake command not found. Please ensure CMake is installed and in your PATH.", file=sys.stderr)
        sys.exit(1)

    # Step 2: Build using CMake
    print("\n--- Building project with CMake ---")
    build_cmd = ["cmake", "--build", str(build_dir)]
    try:
        subprocess.run(build_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: CMake build failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)

    # Step 3: Run the compiled executable
    if not binary_path.exists():
        print(f"Error: Executable not found at {binary_path}", file=sys.stderr)
        sys.exit(1)

    print("\n--- Running coherent_thomson_solver ---")
    run_cmd = [str(binary_path)]
    if args.config_file:
        # Resolve config file path (relative to current working directory)
        config_path = Path(args.config_file).resolve()
        if not config_path.exists():
            print(f"Warning: Configuration file '{args.config_file}' not found.", file=sys.stderr)
        run_cmd.append(str(config_path))
        print(f"Passing configuration file: {config_path}")

    try:
        result = subprocess.run(run_cmd)
        sys.exit(result.returncode)
    except Exception as e:
        print(f"Error executing solver: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
