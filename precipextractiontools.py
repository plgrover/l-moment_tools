import pandas as pd

'''
Extracts the rainfall from the database
'''
def getPrecipitationRecord(stationId, conn):
    sql = 'SELECT reading_dt,accumulation_mm from precipitation.station_precipitation \
        where climate_stations_id={0} AND duration_hrs=24;'.format(stationId)

    dailyprecipDF = pd.read_sql(sql, conn)
    dailyprecipDF = dailyprecipDF.set_index('reading_dt')
    return dailyprecipDF


def findNearbyStation(lat, long, conn, range_deg=0.25):
    sql = 'SELECT * from precipitation.climate_stations WHERE latitude between {0} and {1} AND longitude between {2} and {3};'
    sql = sql.format(lat - range_deg, lat + range_deg, long - range_deg, long + range_deg)
    nearbyStations = pd.read_sql(sql, conn)
    return nearbyStations


'''
Identifies years with bad data.
'''
def getBadYearsDaily(dailyprecipDF, percent_missing=0.1):
    daysPerYear = 365
    countTs = dailyprecipDF.groupby(dailyprecipDF.index.year).count()
    retval = countTs[countTs.accumulation_mm < (1. - percent_missing) * daysPerYear]
    return retval


'''
Generates the annual maximum series for a specified duration.
'''
def getAnnualMaximumSeriesforDuration(dailyprecipDF, badYears=None, duration='7d'):
    durationdf = dailyprecipDF.rolling(window=duration, min_periods=1).sum()
    aMaxDurationDf = durationdf.groupby(durationdf.index.year).max()
    if badYears is not None:
        aMaxDurationDf = aMaxDurationDf.drop(labels=badYears.index)
    return aMaxDurationDf