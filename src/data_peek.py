import pandas as pd

def peekTable(df, size=5):
    """
    Prints the first n (default=5) rows and columns of a pandas DataFrame.

    Args:
        df: pandas DataFrame to peek at.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")

    num_rows = min(size, len(df))
    num_cols = min(size, len(df.columns))

    if num_rows == 0:
        print("DataFrame is empty.")
        return

    peek = df.iloc[:num_rows, :num_cols]
    print(peek)
