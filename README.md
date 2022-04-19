## CH3OrderParams
### Python script to calculate methyl order parameters in proteins

*Requieres*:
- Pandas 
- argparse

Takes as input amber NetCDF trajectories and topology to and returns a CSV file with order parameters of methyl groups calculated as stated in https://www.cell.com/action/showPdf?pii=S0006-3495%2811%2900783-1

- -t --topology String containing topology for MD trajectory
- -x --trajectory Python list of trayectory names as list of strings
- -o --output Output name for CSV fileCancel changes
- -r --residues default:'ILVAM' "Methyl aminoacids"
