import argparse

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Process a CSV file and optional list of strings.')

    # Add the required argument for the CSV file
    parser.add_argument('csv_file', type=str, help='Path to the input CSV file')

    # Add the optional argument for a list of strings
    parser.add_argument('--strings',"-s", type=str, nargs='*', help='Optional list of strings')
    parser.add_argument("--only_prot", action="store_true")

    # Parse the arguments
    args = parser.parse_args()

    # Print the arguments for demonstration
    print(f'CSV File: {args.csv_file}')
    if args.strings:
        print(f'List of Strings: {args.strings}')
    else:
        print('No additional strings provided.')
    print({args.strings})

if __name__ == '__main__':
    main()
