"""common custom Python functions for use in other scripts"""

#transpose and set first row as df header
def flip_df(df):
	df=df.T.reset_index()
	df.columns = df.iloc[0]
	df=df.iloc[1:]
	return df

#make last column of dataframe first column
def swap_cols(df):
	cols = df.columns.tolist()
	cols = cols[-1:] + cols[:-1]
	df = df[cols]
	return df

