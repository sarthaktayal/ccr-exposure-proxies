"""Regenerate every figure in ./figures by running both studies."""
import runpy, os
here = os.path.dirname(__file__)
for exp in ("part1_linear.py", "part2_nonlinear.py"):
    print(f"\n=== running {exp} ===")
    runpy.run_path(os.path.join(here, "experiments", exp), run_name="__main__")
print("\nAll figures regenerated in", os.path.join(here, "figures"))
