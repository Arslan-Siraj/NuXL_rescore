from .cli import build_parser, dispatch

def main():
    parser = build_parser()
    args = parser.parse_args() # args come from command line
    dispatch(args)
    