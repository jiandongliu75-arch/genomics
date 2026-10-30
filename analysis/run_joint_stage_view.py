"""Run one independent spatial view in an isolated staging directory."""
import argparse
from pathlib import Path
import joint_signaling_common as common

def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=common.STAGES,required=True)
 p.add_argument('--view',choices=['core','native'],required=True);p.add_argument('--staging',required=True)
 a=p.parse_args();out=Path(a.staging).resolve()
 assert out.is_relative_to(common.OUT) and out!=common.OUT
 common.OUT=out
 import joint_signaling_spatial as spatial
 spatial.prepare(a.stage,[a.view])
if __name__=='__main__':main()
