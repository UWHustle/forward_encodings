CREATE TABLE part
(
    p_partkey   UINTEGER,
    p_name      VARCHAR,
    p_mfgr      VARCHAR,
    p_category  VARCHAR,
    p_brand1    VARCHAR,
    p_color     VARCHAR,
    p_type      VARCHAR,
    p_size      UTINYINT,
    p_container VARCHAR,
    PRIMARY KEY (p_partkey)
);

CREATE TABLE supplier
(
    s_suppkey UINTEGER,
    s_name    VARCHAR,
    s_address VARCHAR,
    s_city    VARCHAR,
    s_nation  VARCHAR,
    s_region  VARCHAR,
    s_phone   VARCHAR,
    PRIMARY KEY (s_suppkey)
);

CREATE TABLE customer
(
    c_custkey    UINTEGER,
    c_name       VARCHAR,
    c_address    VARCHAR,
    c_city       VARCHAR,
    c_nation     VARCHAR,
    c_region     VARCHAR,
    c_phone      VARCHAR,
    c_mktsegment VARCHAR,
    PRIMARY KEY (c_custkey)
);

CREATE TABLE date
(
    d_datekey          UINTEGER,
    d_date             VARCHAR,
    d_dayofweek        VARCHAR,
    d_month            VARCHAR,
    d_year             USMALLINT,
    d_yearmonthnum     UINTEGER,
    d_yearmonth        VARCHAR,
    d_daynuminweek     UTINYINT,
    d_daynuminmonth    UTINYINT,
    d_daynuminyear     USMALLINT,
    d_monthnuminyear   UTINYINT,
    d_weeknuminyear    UTINYINT,
    d_sellingseason    VARCHAR,
    d_lastdayinweekfl  BOOLEAN,
    d_lastdayinmonthfl BOOLEAN,
    d_holidayfl        BOOLEAN,
    d_weekdayfl        BOOLEAN,
    PRIMARY KEY (d_datekey)
);

CREATE TABLE lineorder
(
    lo_orderkey      UINTEGER,
    lo_linenumber    UTINYINT,
    lo_custkey       UINTEGER,
    lo_partkey       UINTEGER,
    lo_suppkey       UINTEGER,
    lo_orderdate     UINTEGER,
    lo_orderpriority VARCHAR,
    lo_shippriority  VARCHAR,
    lo_quantity      UTINYINT,
    lo_extendedprice UINTEGER,
    lo_ordtotalprice UINTEGER,
    lo_discount      UTINYINT,
    lo_revenue       UINTEGER,
    lo_supplycost    UINTEGER,
    lo_tax           UTINYINT,
    lo_commitdate    UINTEGER,
    lo_shipmode      VARCHAR,
    PRIMARY KEY (lo_orderkey, lo_linenumber)
);

COPY part FROM '?1ssb-dbgen/build/?2part.tbl' ( DELIMITER '|' );
COPY supplier FROM '?1ssb-dbgen/build/?2supplier.tbl' ( DELIMITER '|' );
COPY customer FROM '?1ssb-dbgen/build/?2customer.tbl' ( DELIMITER '|' );
COPY date FROM '?1ssb-dbgen/build/?2date.tbl' ( DELIMITER '|' );
COPY lineorder FROM '?1ssb-dbgen/build/?2lineorder.tbl' ( DELIMITER '|' );

COPY (SELECT lo_quantity,
        lo_extendedprice,
        lo_discount,
        lo_revenue,
        lo_supplycost,
        p_mfgr,
        p_category,
        p_brand1,
        s_city,
        s_nation,
        s_region,
        c_city,
        c_nation,
        c_region,
        d_year,
        d_yearmonthnum,
        d_yearmonth,
        d_weeknuminyear
      FROM lineorder, part, supplier, customer, date
      WHERE lo_partkey = p_partkey
        AND lo_suppkey = s_suppkey
        AND lo_custkey = c_custkey
        AND lo_orderdate = d_datekey)
TO '?1?2ssb_wide_table.tbl' WITH (DELIMITER '|');
