"""Portable report layout, with only the21-gene maturation outcome."""
from pathlib import Path
import write_morphology_nonCM_report as writer
if __name__=='__main__':
 writer.OUT=Path(__file__).resolve().parents[1]/'results/CM21_primary_spatial'
 writer.STEM='CM21_maturation_and_neighbors'
 writer.FOOTER_LABEL='POSTNATAL HEART | 21-GENE CM MATURATION'
 writer.METHODS_HEADING='## Methods and scope'
 writer.main()
