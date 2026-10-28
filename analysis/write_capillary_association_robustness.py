"""Reuse the validated document layout for the separate P7 follow-up bundle."""
from pathlib import Path
import write_capillary_maturation_retest as writer
if __name__=='__main__':
 writer.OUT=Path(__file__).resolve().parents[1]/'results/capillary_association_robustness'
 writer.main()
