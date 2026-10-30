"""Word/PDF from the same narrative and integrated figure plates."""
from pathlib import Path
import write_morphology_nonCM_report as writer
if __name__=='__main__':
 writer.OUT=Path(__file__).resolve().parents[1]/'results/CM21_paperstyle_extension'
 writer.STEM='CM21_paperstyle_vascular_stromal_report'
 writer.FOOTER_LABEL='POSTNATAL HEART | MEASURED 21-GENE MATURATION'
 writer.METHODS_HEADING='## Methods and scope'
 writer.main()
