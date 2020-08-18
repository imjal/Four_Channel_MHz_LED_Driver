EESchema Schematic File Version 4
EELAYER 30 0
EELAYER END
$Descr A2 23386 16535
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L Custom_parts:DAC084S085 U?
U 1 1 5F33B1A2
P 2750 5650
F 0 "U?" H 3400 5013 60  0000 C CNN
F 1 "DAC084S085" H 3400 5119 60  0000 C CNN
F 2 "Custom Footprints:DAC084S085CIMM" H 3450 5850 60  0001 C CNN
F 3 "http://www.ti.com/general/docs/suppproductinfo.tsp?distId=10&gotoUrl=http%3A%2F%2Fwww.ti.com%2Flit%2Fgpn%2Fdac084s085" H 3400 5119 60  0001 C CNN
F 4 "Texas Instruments" H 2750 5650 50  0001 C CNN "Manufacturer"
F 5 "DAC084S085CIMM/NOPB" H 2750 5650 50  0001 C CNN "Part #"
	1    2750 5650
	-1   0    0    1   
$EndComp
$Comp
L Custom_parts:MCP4351-502E_ST U?
U 1 1 5F3407BD
P 3950 5000
F 0 "U?" H 4550 5257 60  0000 C CNN
F 1 "MCP4351-502E_ST" H 4550 5151 60  0000 C CNN
F 2 "Custom Footprints:MCP4351-502E" H 4600 5250 60  0001 C CNN
F 3 "http://www.microchip.com/mymicrochip/filehandler.aspx?ddocname=en547555" H 4550 5151 60  0001 C CNN
F 4 "Microchip Technology" H 3950 5000 50  0001 C CNN "Manufacturer"
F 5 "MCP4351-502E/ST" H 3950 5000 50  0001 C CNN "Part #"
	1    3950 5000
	1    0    0    -1  
$EndComp
Wire Wire Line
	5150 5900 5150 6000
Wire Wire Line
	5150 4700 5150 5000
Wire Wire Line
	2750 5250 2900 5250
Wire Wire Line
	2950 5350 2750 5350
$Comp
L power:GND #PWR?
U 1 1 5F34784A
P 3950 5200
F 0 "#PWR?" H 3950 4950 50  0001 C CNN
F 1 "GND" V 3955 5072 50  0000 R CNN
F 2 "" H 3950 5200 50  0001 C CNN
F 3 "" H 3950 5200 50  0001 C CNN
	1    3950 5200
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F34877E
P 3950 5700
F 0 "#PWR?" H 3950 5450 50  0001 C CNN
F 1 "GND" V 3850 5700 50  0000 R CNN
F 2 "" H 3950 5700 50  0001 C CNN
F 3 "" H 3950 5700 50  0001 C CNN
	1    3950 5700
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F348E43
P 5150 5700
F 0 "#PWR?" H 5150 5450 50  0001 C CNN
F 1 "GND" V 5155 5572 50  0000 R CNN
F 2 "" H 5150 5700 50  0001 C CNN
F 3 "" H 5150 5700 50  0001 C CNN
	1    5150 5700
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F349663
P 5150 5200
F 0 "#PWR?" H 5150 4950 50  0001 C CNN
F 1 "GND" V 5155 5072 50  0000 R CNN
F 2 "" H 5150 5200 50  0001 C CNN
F 3 "" H 5150 5200 50  0001 C CNN
	1    5150 5200
	0    -1   -1   0   
$EndComp
Wire Wire Line
	3950 5600 3950 5700
Connection ~ 3950 5700
NoConn ~ 5150 5400
NoConn ~ 5150 5500
Text Label 5150 5300 0    50   ~ 0
5V-analog
$Comp
L power:GND #PWR?
U 1 1 5F34BA64
P 1450 5250
F 0 "#PWR?" H 1450 5000 50  0001 C CNN
F 1 "GND" V 1455 5122 50  0000 R CNN
F 2 "" H 1450 5250 50  0001 C CNN
F 3 "" H 1450 5250 50  0001 C CNN
	1    1450 5250
	0    1    1    0   
$EndComp
Text Label 1450 5350 2    50   ~ 0
5V-analog
$Comp
L power:GND #PWR?
U 1 1 5F3504EC
P 2750 5850
F 0 "#PWR?" H 2750 5600 50  0001 C CNN
F 1 "GND" V 2755 5722 50  0000 R CNN
F 2 "" H 2750 5850 50  0001 C CNN
F 3 "" H 2750 5850 50  0001 C CNN
	1    2750 5850
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F352F78
P 5250 5500
F 0 "#PWR?" H 5250 5250 50  0001 C CNN
F 1 "GND" H 5450 5400 50  0000 R CNN
F 2 "" H 5250 5500 50  0001 C CNN
F 3 "" H 5250 5500 50  0001 C CNN
	1    5250 5500
	1    0    0    -1  
$EndComp
Wire Wire Line
	5250 5300 5150 5300
$Comp
L Custom_parts:S18V20F12_12V_DC U?
U 1 1 5F357FC0
P 2700 1650
F 0 "U?" H 2675 1785 50  0000 C CNN
F 1 "S18V20F12_12V_DC" H 2675 1694 50  0000 C CNN
F 2 "" H 2700 1650 50  0001 C CNN
F 3 "https://www.pololu.com/product-info-merged/2577" H 2700 1650 50  0001 C CNN
F 4 "Pololu Corporation" H 2700 1650 50  0001 C CNN "Manufacturer"
F 5 "2577" H 2700 1650 50  0001 C CNN "Part #"
	1    2700 1650
	1    0    0    -1  
$EndComp
Wire Wire Line
	2200 1750 2200 1800
Wire Wire Line
	2200 1950 2200 2000
Wire Wire Line
	3150 1950 3150 2000
Text Label 2100 2000 2    50   ~ 0
Vin
Connection ~ 2200 2000
Wire Wire Line
	2200 2000 2200 2050
Wire Wire Line
	3150 2000 3250 2000
Connection ~ 3150 2000
Wire Wire Line
	3150 2000 3150 2050
$Comp
L power:GND #PWR?
U 1 1 5F35A16A
P 2200 1800
F 0 "#PWR?" H 2200 1550 50  0001 C CNN
F 1 "GND" V 2205 1672 50  0000 R CNN
F 2 "" H 2200 1800 50  0001 C CNN
F 3 "" H 2200 1800 50  0001 C CNN
	1    2200 1800
	0    1    1    0   
$EndComp
Connection ~ 2200 1800
Wire Wire Line
	2200 1800 2200 1850
Text Label 3250 2000 0    50   ~ 0
12V
$Comp
L Custom_parts:DPBW06F-05 U?
U 1 1 5F35D26E
P 4550 1650
F 0 "U?" H 4600 1785 50  0000 C CNN
F 1 "DPBW06F-05" H 4600 1694 50  0000 C CNN
F 2 "Custom Footprints:DPBW06F-05" H 4550 1650 50  0001 C CNN
F 3 "https://media.digikey.com/pdf/Data%20Sheets/Mean%20Well%20PDF's/SPBW06,DPBW06_Ds.pdf" H 4550 1650 50  0001 C CNN
F 4 "MEAN WELL USA Inc." H 4550 1650 50  0001 C CNN "Manufacturer"
F 5 "DPBW06F-05" H 4550 1650 50  0001 C CNN "Part #"
	1    4550 1650
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F35E7FC
P 3250 1850
F 0 "C?" H 3342 1896 50  0000 L CNN
F 1 "47uF" H 3342 1805 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 3250 1850 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 3250 1850 50  0001 C CNN
F 4 "Murata Electronics" H 3250 1850 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 3250 1850 50  0001 C CNN "Part #"
	1    3250 1850
	1    0    0    -1  
$EndComp
NoConn ~ 4100 1950
$Comp
L power:GND #PWR?
U 1 1 5F35F7A3
P 3250 1750
F 0 "#PWR?" H 3250 1500 50  0001 C CNN
F 1 "GND" H 3150 1650 50  0000 R CNN
F 2 "" H 3250 1750 50  0001 C CNN
F 3 "" H 3250 1750 50  0001 C CNN
	1    3250 1750
	-1   0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F3602CD
P 4100 1750
F 0 "#PWR?" H 4100 1500 50  0001 C CNN
F 1 "GND" V 4105 1622 50  0000 R CNN
F 2 "" H 4100 1750 50  0001 C CNN
F 3 "" H 4100 1750 50  0001 C CNN
	1    4100 1750
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F36457B
P 5550 1850
F 0 "#PWR?" H 5550 1600 50  0001 C CNN
F 1 "GND" V 5555 1722 50  0000 R CNN
F 2 "" H 5550 1850 50  0001 C CNN
F 3 "" H 5550 1850 50  0001 C CNN
	1    5550 1850
	0    -1   -1   0   
$EndComp
$Comp
L pspice:INDUCTOR L?
U 1 1 5F3651AC
P 3650 2000
F 0 "L?" H 3650 2215 50  0000 C CNN
F 1 "10uH" H 3650 2124 50  0000 C CNN
F 2 "Inductor_SMD:L_1210_3225Metric" H 3650 2000 50  0001 C CNN
F 3 "https://product.tdk.com/info/en/catalog/datasheets/inductor_automotive_power_tfm322512alma_en.pdf?ref_disty=digikey" H 3650 2000 50  0001 C CNN
F 4 "TDK Corporation" H 3650 2000 50  0001 C CNN "Manufacturer"
F 5 "TFM322512ALMA100MTAA" H 3650 2000 50  0001 C CNN "Part #"
	1    3650 2000
	1    0    0    -1  
$EndComp
Wire Wire Line
	3400 2000 3250 2000
Connection ~ 3250 2000
Wire Wire Line
	3900 2000 3900 1850
Wire Wire Line
	3900 1850 4100 1850
$Comp
L Device:C_Small C?
U 1 1 5F36881C
P 4600 2200
F 0 "C?" V 4500 2200 50  0000 C CNN
F 1 "1500pF" V 4400 2200 50  0000 C CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 4600 2200 50  0001 C CNN
F 3 "https://katalog.we-online.de/pbs/datasheet/885342206004.pdf" H 4600 2200 50  0001 C CNN
F 4 "Würth Elektronik" H 4600 2200 50  0001 C CNN "Manufacturer"
F 5 "885342206004" H 4600 2200 50  0001 C CNN "Part #"
	1    4600 2200
	0    1    -1   0   
$EndComp
Wire Wire Line
	4500 2200 3900 2200
Wire Wire Line
	3900 2200 3900 2000
Connection ~ 3900 2000
Wire Wire Line
	4700 2200 5100 2200
$Comp
L Device:C_Small C?
U 1 1 5F36C3BB
P 4600 1400
F 0 "C?" V 4700 1400 50  0000 C CNN
F 1 "1500pF" V 4650 1600 50  0000 C CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 4600 1400 50  0001 C CNN
F 3 "https://katalog.we-online.de/pbs/datasheet/885342206004.pdf" H 4600 1400 50  0001 C CNN
F 4 "Würth Elektronik" H 4600 1400 50  0001 C CNN "Manufacturer"
F 5 "885342206004" H 4600 1400 50  0001 C CNN "Part #"
	1    4600 1400
	0    1    -1   0   
$EndComp
Wire Wire Line
	4100 1750 4100 1400
Wire Wire Line
	4100 1400 4500 1400
Connection ~ 4100 1750
Wire Wire Line
	5100 1750 5100 1650
Wire Wire Line
	5100 1400 4700 1400
$Comp
L Device:C_Small C?
U 1 1 5F370F72
P 5300 1750
F 0 "C?" H 5392 1796 50  0000 L CNN
F 1 "47uF" H 5392 1705 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 5300 1750 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 5300 1750 50  0001 C CNN
F 4 "Murata Electronics" H 5300 1750 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 5300 1750 50  0001 C CNN "Part #"
	1    5300 1750
	1    0    0    -1  
$EndComp
Connection ~ 5100 1650
Wire Wire Line
	5100 1650 5100 1400
Text Label 5150 2050 0    50   ~ 0
5V
Text Label 5100 1650 0    50   ~ 0
-5V
Connection ~ 5300 1850
Wire Wire Line
	5100 1850 5300 1850
Wire Wire Line
	5100 1650 5300 1650
Wire Wire Line
	3250 1950 3250 2000
Wire Wire Line
	3150 1750 3150 1850
Wire Wire Line
	3150 1750 3250 1750
Connection ~ 3150 1750
Connection ~ 3250 1750
$Comp
L Device:R_Small R?
U 1 1 5F379C80
P 3250 2100
F 0 "R?" H 3309 2146 50  0000 L CNN
F 1 "100" H 3309 2055 50  0000 L CNN
F 2 "Resistor_SMD:R_0805_2012Metric" H 3250 2100 50  0001 C CNN
F 3 "https://www.seielect.com/catalog/sei-rncp.pdf" H 3250 2100 50  0001 C CNN
F 4 "Stackpole Electronics Inc" H 3250 2100 50  0001 C CNN "Manufacturer"
F 5 "RNCP0805FTD100R" H 3250 2100 50  0001 C CNN "Part #"
	1    3250 2100
	1    0    0    -1  
$EndComp
$Comp
L Custom_parts:ADP7118AUJZ-5.0-R7 U?
U 1 1 5F380EEB
P 3550 2700
F 0 "U?" H 3975 2907 60  0000 C CNN
F 1 "ADP7118AUJZ-5.0-R7" H 3975 2801 60  0000 C CNN
F 2 "Custom Footprints:ADP7118AUJZ-R7" H 4650 2940 60  0001 C CNN
F 3 "https://www.analog.com/media/en/technical-documentation/data-sheets/ADP7118.pdf" H 3975 2801 60  0001 C CNN
F 4 "Analog Devices Inc." H 3550 2700 50  0001 C CNN "Manufacturer"
F 5 "ADP7118AUJZ-5.0-R7" H 3550 2700 50  0001 C CNN "Part #"
	1    3550 2700
	1    0    0    -1  
$EndComp
Connection ~ 3250 2750
Wire Wire Line
	4700 2950 4700 2750
Text Label 4800 2750 0    50   ~ 0
5V-analog
Wire Wire Line
	3250 2200 3250 2350
$Comp
L Device:C_Small C?
U 1 1 5F388B05
P 4900 2850
F 0 "C?" H 4992 2896 50  0000 L CNN
F 1 "2.2uF" H 4992 2805 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 4900 2850 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 4900 2850 50  0001 C CNN
F 4 "Taiyo Yuden" H 4900 2850 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 4900 2850 50  0001 C CNN "Part #"
	1    4900 2850
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F38AF36
P 4900 2950
F 0 "#PWR?" H 4900 2700 50  0001 C CNN
F 1 "GND" H 5100 2900 50  0000 R CNN
F 2 "" H 4900 2950 50  0001 C CNN
F 3 "" H 4900 2950 50  0001 C CNN
	1    4900 2950
	1    0    0    -1  
$EndComp
Wire Wire Line
	4900 2750 4700 2750
Connection ~ 4700 2750
$Comp
L Device:C_Small C?
U 1 1 5F38C930
P 3150 2750
F 0 "C?" V 3350 2700 50  0000 L CNN
F 1 "2.2uF" V 3250 2650 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 3150 2750 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 3150 2750 50  0001 C CNN
F 4 "Taiyo Yuden" H 3150 2750 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 3150 2750 50  0001 C CNN "Part #"
	1    3150 2750
	0    -1   -1   0   
$EndComp
Connection ~ 3250 2350
$Comp
L power:GND #PWR?
U 1 1 5F39026E
P 3050 2350
F 0 "#PWR?" H 3050 2100 50  0001 C CNN
F 1 "GND" V 3050 2200 50  0000 R CNN
F 2 "" H 3050 2350 50  0001 C CNN
F 3 "" H 3050 2350 50  0001 C CNN
	1    3050 2350
	0    1    1    0   
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F390986
P 5250 5400
F 0 "C?" H 5342 5446 50  0000 L CNN
F 1 "2.2uF" H 5342 5355 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 5250 5400 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 5250 5400 50  0001 C CNN
F 4 "Taiyo Yuden" H 5250 5400 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 5250 5400 50  0001 C CNN "Part #"
	1    5250 5400
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F3918D7
P 2750 5750
F 0 "C?" H 2550 5750 50  0000 L CNN
F 1 "2.2uF" H 2500 5650 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 2750 5750 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 2750 5750 50  0001 C CNN
F 4 "Taiyo Yuden" H 2750 5750 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 2750 5750 50  0001 C CNN "Part #"
	1    2750 5750
	1    0    0    -1  
$EndComp
Text Label 2750 5650 0    50   ~ 0
5V-analog
Wire Wire Line
	2750 5450 3250 5450
Wire Wire Line
	2750 5550 3200 5550
Text Notes 2800 3100 0    59   ~ 0
LDO - 12V to clean 5V for analog circuits
Text Notes 3800 1250 0    59   ~ 0
ISOLATED - 12V to split +/- 5V
Text Notes 1800 1500 0    59   ~ 0
SEPIC - Vin (3V - 30V) to 12V DC
Wire Wire Line
	3250 2350 3250 2450
$Comp
L Jumper:SolderJumper_2_Open JP?
U 1 1 5F3A0F69
P 2200 2350
F 0 "JP?" H 2200 2555 50  0000 C CNN
F 1 "SolderJumper_2_Open" H 2250 2450 50  0000 C CNN
F 2 "Jumper:SolderJumper-2_P1.3mm_Open_Pad1.0x1.5mm" H 2200 2350 50  0001 C CNN
F 3 "~" H 2200 2350 50  0001 C CNN
	1    2200 2350
	1    0    0    -1  
$EndComp
Text Label 2350 2350 0    50   ~ 0
12V
$Comp
L Device:C_Small C?
U 1 1 5F3A29AA
P 3150 2350
F 0 "C?" V 2921 2350 50  0000 C CNN
F 1 "47uF" V 3012 2350 50  0000 C CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 3150 2350 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 3150 2350 50  0001 C CNN
F 4 "Murata Electronics" H 3150 2350 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 3150 2350 50  0001 C CNN "Part #"
	1    3150 2350
	0    1    1    0   
$EndComp
Connection ~ 3250 2450
Wire Wire Line
	3250 2450 3250 2750
$Comp
L power:GND #PWR?
U 1 1 5F3A55D8
P 3050 2750
F 0 "#PWR?" H 3050 2500 50  0001 C CNN
F 1 "GND" V 3050 2600 50  0000 R CNN
F 2 "" H 3050 2750 50  0001 C CNN
F 3 "" H 3050 2750 50  0001 C CNN
	1    3050 2750
	0    1    1    0   
$EndComp
Wire Wire Line
	3050 2750 3050 2850
Wire Wire Line
	3050 2850 3250 2850
Connection ~ 3050 2750
Wire Wire Line
	2700 2950 2700 2450
Wire Wire Line
	2700 2950 3250 2950
Wire Wire Line
	2700 2450 3250 2450
Text Notes 1750 2500 0    59   ~ 0
12V Supply Bypass
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3BA201
P 7150 4000
F 0 "U?" H 7325 4165 50  0000 C CNN
F 1 "TMUX1204DGSR" H 7325 4074 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 7150 5000 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 7150 4000 50  0001 C CNN
F 4 "Texas Instruments" H 7150 4000 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 7150 4000 50  0001 C CNN "Part #"
	1    7150 4000
	1    0    0    -1  
$EndComp
Wire Wire Line
	6950 4450 6950 4550
Wire Wire Line
	6950 4550 7700 4550
Wire Wire Line
	7700 4550 7700 4450
Text Label 6950 4050 2    50   ~ 0
interline_1
Wire Wire Line
	7700 4150 7750 4150
Wire Wire Line
	7750 4150 7750 4350
Wire Wire Line
	7750 4350 7700 4350
Text Label 7700 4050 0    50   ~ 0
analog_select_1
$Comp
L power:GND #PWR?
U 1 1 5F3BE358
P 6950 4250
F 0 "#PWR?" H 6950 4000 50  0001 C CNN
F 1 "GND" V 6955 4122 50  0000 R CNN
F 2 "" H 6950 4250 50  0001 C CNN
F 3 "" H 6950 4250 50  0001 C CNN
	1    6950 4250
	0    1    1    0   
$EndComp
Text Label 5150 5800 0    50   ~ 0
internal_analog_1
Text Label 5150 5100 0    50   ~ 0
internal_analog_3
Text Label 3950 5100 2    50   ~ 0
internal_analog_4
Text Label 3950 5800 2    50   ~ 0
internal_analog_2
Wire Wire Line
	2900 5000 2900 5250
Wire Wire Line
	2900 5000 3950 5000
Wire Wire Line
	2950 4700 2950 5350
Wire Wire Line
	2950 4700 5150 4700
Wire Wire Line
	3250 5900 3250 5450
Wire Wire Line
	3250 5900 3950 5900
Wire Wire Line
	3200 5550 3200 6000
Wire Wire Line
	3200 6000 5150 6000
Text Label 6950 4150 2    50   ~ 0
internal_analog_1
Text Label 6950 4350 2    50   ~ 0
external_analog_1
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3DA347
P 7150 4850
F 0 "U?" H 7325 5015 50  0000 C CNN
F 1 "TMUX1204DGSR" H 7325 4924 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 7150 5850 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 7150 4850 50  0001 C CNN
F 4 "Texas Instruments" H 7150 4850 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 7150 4850 50  0001 C CNN "Part #"
	1    7150 4850
	1    0    0    -1  
$EndComp
Wire Wire Line
	6950 5300 6950 5400
Wire Wire Line
	6950 5400 7700 5400
Wire Wire Line
	7700 5400 7700 5300
Text Label 6950 4900 2    50   ~ 0
interline_2
Wire Wire Line
	7700 5000 7750 5000
Wire Wire Line
	7750 5000 7750 5200
Wire Wire Line
	7750 5200 7700 5200
Text Label 7700 4900 0    50   ~ 0
analog_select_2
$Comp
L power:GND #PWR?
U 1 1 5F3DA355
P 6950 5100
F 0 "#PWR?" H 6950 4850 50  0001 C CNN
F 1 "GND" V 6955 4972 50  0000 R CNN
F 2 "" H 6950 5100 50  0001 C CNN
F 3 "" H 6950 5100 50  0001 C CNN
	1    6950 5100
	0    1    1    0   
$EndComp
Text Label 6950 5000 2    50   ~ 0
internal_analog_2
Text Label 6950 5200 2    50   ~ 0
external_analog_2
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3DC268
P 7200 5700
F 0 "U?" H 7375 5865 50  0000 C CNN
F 1 "TMUX1204DGSR" H 7375 5774 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 7200 6700 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 7200 5700 50  0001 C CNN
F 4 "Texas Instruments" H 7200 5700 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 7200 5700 50  0001 C CNN "Part #"
	1    7200 5700
	1    0    0    -1  
$EndComp
Wire Wire Line
	7000 6150 7000 6250
Wire Wire Line
	7000 6250 7750 6250
Wire Wire Line
	7750 6250 7750 6150
Text Label 7000 5750 2    50   ~ 0
interline_3
Wire Wire Line
	7750 5850 7800 5850
Wire Wire Line
	7800 5850 7800 6050
Wire Wire Line
	7800 6050 7750 6050
Text Label 7750 5750 0    50   ~ 0
analog_select_3
$Comp
L power:GND #PWR?
U 1 1 5F3DC276
P 7000 5950
F 0 "#PWR?" H 7000 5700 50  0001 C CNN
F 1 "GND" V 7005 5822 50  0000 R CNN
F 2 "" H 7000 5950 50  0001 C CNN
F 3 "" H 7000 5950 50  0001 C CNN
	1    7000 5950
	0    1    1    0   
$EndComp
Text Label 7000 5850 2    50   ~ 0
internal_analog_3
Text Label 7000 6050 2    50   ~ 0
external_analog_3
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3DEF78
P 7200 6500
F 0 "U?" H 7375 6665 50  0000 C CNN
F 1 "TMUX1204DGSR" H 7375 6574 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 7200 7500 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 7200 6500 50  0001 C CNN
F 4 "Texas Instruments" H 7200 6500 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 7200 6500 50  0001 C CNN "Part #"
	1    7200 6500
	1    0    0    -1  
$EndComp
Wire Wire Line
	7000 6950 7000 7050
Wire Wire Line
	7000 7050 7750 7050
Wire Wire Line
	7750 7050 7750 6950
Text Label 7000 6550 2    50   ~ 0
interline_4
Wire Wire Line
	7750 6650 7800 6650
Wire Wire Line
	7800 6650 7800 6850
Wire Wire Line
	7800 6850 7750 6850
Text Label 7750 6550 0    50   ~ 0
analog_select_4
$Comp
L power:GND #PWR?
U 1 1 5F3DEF86
P 7000 6750
F 0 "#PWR?" H 7000 6500 50  0001 C CNN
F 1 "GND" V 7005 6622 50  0000 R CNN
F 2 "" H 7000 6750 50  0001 C CNN
F 3 "" H 7000 6750 50  0001 C CNN
	1    7000 6750
	0    1    1    0   
$EndComp
Text Label 7000 6650 2    50   ~ 0
internal_analog_4
Text Label 7000 6850 2    50   ~ 0
external_analog_4
Text Notes 2000 4550 0    59   ~ 0
Internal 4-channel AWG with current limit
$Comp
L Connector:Conn_Coaxial_x2 J?
U 1 1 5F3FAF61
P 1350 7650
F 0 "J?" H 1300 8000 50  0000 L CNN
F 1 "Conn_Coaxial_x2" H 1050 7900 50  0000 L CNN
F 2 "Custom Footprints:031-6575_2x_BNC" H 1350 7550 50  0001 C CNN
F 3 "https://www.amphenolrf.com/library/download/link/link_id/436470/parent/031-6575/" H 1350 7550 50  0001 C CNN
F 4 "Amphenol RF" H 1350 7650 50  0001 C CNN "Manufacturer"
F 5 "031-6575" H 1350 7650 50  0001 C CNN "Part #"
	1    1350 7650
	-1   0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F3FBEF1
P 1350 7050
F 0 "#PWR?" H 1350 6800 50  0001 C CNN
F 1 "GND" H 1450 6900 50  0000 R CNN
F 2 "" H 1350 7050 50  0001 C CNN
F 3 "" H 1350 7050 50  0001 C CNN
	1    1350 7050
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F3FC8FB
P 1350 7950
F 0 "#PWR?" H 1350 7700 50  0001 C CNN
F 1 "GND" H 1450 7800 50  0000 R CNN
F 2 "" H 1350 7950 50  0001 C CNN
F 3 "" H 1350 7950 50  0001 C CNN
	1    1350 7950
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_Coaxial_x2 J?
U 1 1 5F40799C
P 1350 6750
F 0 "J?" H 1300 7100 50  0000 L CNN
F 1 "Conn_Coaxial_x2" H 1050 7000 50  0000 L CNN
F 2 "Custom Footprints:031-6575_2x_BNC" H 1350 6650 50  0001 C CNN
F 3 "https://www.amphenolrf.com/library/download/link/link_id/436470/parent/031-6575/" H 1350 6650 50  0001 C CNN
F 4 "Amphenol RF" H 1350 6750 50  0001 C CNN "Manufacturer"
F 5 "031-6575" H 1350 6750 50  0001 C CNN "Part #"
	1    1350 6750
	-1   0    0    -1  
$EndComp
Text Label 2600 7100 2    50   ~ 0
external_analog_1
Text Label 2400 7850 0    50   ~ 0
external_analog_2
Text Label 2400 7950 0    50   ~ 0
external_analog_3
Text Label 2400 8050 0    50   ~ 0
external_analog_4
$Comp
L Switch:SW_DIP_x03 SW?
U 1 1 5F411870
P 5200 7400
F 0 "SW?" H 5200 7867 50  0000 C CNN
F 1 "SW_DIP_x03" H 5200 7776 50  0000 C CNN
F 2 "Custom Footprints:SW_DS04-254-2-03BK-SMT" H 5200 7400 50  0001 C CNN
F 3 "https://www.cuidevices.com/api/videos/videoplayer/smallplayer/ds04-254-smt.pdf" H 5200 7400 50  0001 C CNN
F 4 "CUI Devices" H 5200 7400 50  0001 C CNN "Manufacturer"
F 5 "DS04-254-2-03BK-SMT" H 5200 7400 50  0001 C CNN "Part #"
	1    5200 7400
	1    0    0    -1  
$EndComp
Connection ~ 5100 2050
Wire Wire Line
	5100 2050 5100 1950
Wire Wire Line
	5100 2200 5100 2050
Wire Wire Line
	5300 1850 5550 1850
Text Label 6000 2150 2    50   ~ 0
-0.25V_analog
$Comp
L power:GND #PWR?
U 1 1 5F3D1FD7
P 6000 2400
F 0 "#PWR?" H 6000 2150 50  0001 C CNN
F 1 "GND" V 6100 2350 50  0000 R CNN
F 2 "" H 6000 2400 50  0001 C CNN
F 3 "" H 6000 2400 50  0001 C CNN
	1    6000 2400
	1    0    0    -1  
$EndComp
$Comp
L Device:R_Small R?
U 1 1 5F3D1763
P 6000 2300
F 0 "R?" H 6059 2346 50  0000 L CNN
F 1 "249" H 6059 2255 50  0000 L CNN
F 2 "Resistor_SMD:R_0603_1608Metric" H 6000 2300 50  0001 C CNN
F 3 "https://www.seielect.com/api/videos/videoplayer/smallplayer/sei-rncp.pdf" H 6000 2300 50  0001 C CNN
F 4 "Stackpole Electronics Inc" H 6000 2300 50  0001 C CNN "Manufacturer"
F 5 "RNCP0603FTD249R" H 6000 2300 50  0001 C CNN "Part #"
	1    6000 2300
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F3D013F
P 6200 1950
F 0 "#PWR?" H 6200 1700 50  0001 C CNN
F 1 "GND" V 6300 1900 50  0000 R CNN
F 2 "" H 6200 1950 50  0001 C CNN
F 3 "" H 6200 1950 50  0001 C CNN
	1    6200 1950
	0    -1   -1   0   
$EndComp
Wire Wire Line
	6000 1950 6000 2200
Connection ~ 6000 1950
Wire Wire Line
	6000 1850 6000 1950
$Comp
L Device:C_Small C?
U 1 1 5F3CD242
P 6100 1950
F 0 "C?" V 6200 1900 50  0000 L CNN
F 1 "47uF" V 6300 1900 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 6100 1950 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 6100 1950 50  0001 C CNN
F 4 "Murata Electronics" H 6100 1950 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 6100 1950 50  0001 C CNN "Part #"
	1    6100 1950
	0    1    1    0   
$EndComp
$Comp
L Device:R_Small R?
U 1 1 5F3CC1B4
P 6000 1750
F 0 "R?" H 6059 1796 50  0000 L CNN
F 1 "4.7k" H 6059 1705 50  0000 L CNN
F 2 "Resistor_SMD:R_0805_2012Metric" H 6000 1750 50  0001 C CNN
F 3 "https://www.seielect.com/catalog/sei-rncp.pdf" H 6000 1750 50  0001 C CNN
F 4 "Yageo" H 6000 1750 50  0001 C CNN "Manufacturer"
F 5 "RT0805FRE074K7L" H 6000 1750 50  0001 C CNN "Part #"
	1    6000 1750
	1    0    0    -1  
$EndComp
Wire Wire Line
	5100 2050 5300 2050
$Comp
L Device:C_Small C?
U 1 1 5F37216A
P 5300 1950
F 0 "C?" H 5392 1996 50  0000 L CNN
F 1 "47uF" H 5392 1905 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 5300 1950 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 5300 1950 50  0001 C CNN
F 4 "Murata Electronics" H 5300 1950 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 5300 1950 50  0001 C CNN "Part #"
	1    5300 1950
	1    0    0    -1  
$EndComp
Wire Wire Line
	5300 1650 6000 1650
Connection ~ 5300 1650
Text Label 7750 4150 0    50   ~ 0
-0.25V_analog
Text Label 7750 5000 0    50   ~ 0
-0.25V_analog
Text Label 7800 5850 0    50   ~ 0
-0.25V_analog
Text Label 7800 6650 0    50   ~ 0
-0.25V_analog
Text Label 7700 4550 0    50   ~ 0
5V-analog
Text Label 7700 5400 0    50   ~ 0
5V-analog
Text Label 7750 6250 0    50   ~ 0
5V-analog
Text Label 7750 7050 0    50   ~ 0
5V-analog
$Comp
L Connector:Conn_Coaxial_x2 J?
U 1 1 5F43D90E
P 1350 9550
F 0 "J?" H 1300 9900 50  0000 L CNN
F 1 "Conn_Coaxial_x2" H 1050 9800 50  0000 L CNN
F 2 "Custom Footprints:031-6575_2x_BNC" H 1350 9450 50  0001 C CNN
F 3 "https://www.amphenolrf.com/library/download/link/link_id/436470/parent/031-6575/" H 1350 9450 50  0001 C CNN
F 4 "Amphenol RF" H 1350 9550 50  0001 C CNN "Manufacturer"
F 5 "031-6575" H 1350 9550 50  0001 C CNN "Part #"
	1    1350 9550
	-1   0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F43D914
P 1350 8950
F 0 "#PWR?" H 1350 8700 50  0001 C CNN
F 1 "GND" H 1450 8800 50  0000 R CNN
F 2 "" H 1350 8950 50  0001 C CNN
F 3 "" H 1350 8950 50  0001 C CNN
	1    1350 8950
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F43D91A
P 1350 9850
F 0 "#PWR?" H 1350 9600 50  0001 C CNN
F 1 "GND" H 1450 9700 50  0000 R CNN
F 2 "" H 1350 9850 50  0001 C CNN
F 3 "" H 1350 9850 50  0001 C CNN
	1    1350 9850
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_Coaxial_x2 J?
U 1 1 5F43D922
P 1350 8650
F 0 "J?" H 1300 9000 50  0000 L CNN
F 1 "Conn_Coaxial_x2" H 1050 8900 50  0000 L CNN
F 2 "Custom Footprints:031-6575_2x_BNC" H 1350 8550 50  0001 C CNN
F 3 "https://www.amphenolrf.com/library/download/link/link_id/436470/parent/031-6575/" H 1350 8550 50  0001 C CNN
F 4 "Amphenol RF" H 1350 8650 50  0001 C CNN "Manufacturer"
F 5 "031-6575" H 1350 8650 50  0001 C CNN "Part #"
	1    1350 8650
	-1   0    0    -1  
$EndComp
Text Label 2400 8150 0    50   ~ 0
A-D_IO_1
Text Label 2400 8250 0    50   ~ 0
A-D_IO_2
Text Label 2400 8350 0    50   ~ 0
A-D_IO_3
Text Label 2400 8450 0    50   ~ 0
A-D_IO_4
$Comp
L Custom_parts:BAT54SDW D?
U 1 1 5F46875C
P 2600 9100
F 0 "D?" H 2775 9447 60  0000 C CNN
F 1 "BAT54SDW" H 2775 9341 60  0000 C CNN
F 2 "Package_TO_SOT_SMD:SOT-363_SC-70-6" H 2800 9300 60  0001 L CNN
F 3 "https://www.diodes.com/assets/Datasheets/ds11005.pdf" H 2800 9400 60  0001 L CNN
F 4 "Diodes Incorporated" H 2600 9100 50  0001 C CNN "Manufacturer"
F 5 "BAT54SDW-7-F" H 2600 9100 50  0001 C CNN "Part #"
	1    2600 9100
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F475E13
P 2400 9000
F 0 "#PWR?" H 2400 8750 50  0001 C CNN
F 1 "GND" V 2300 8950 50  0000 R CNN
F 2 "" H 2400 9000 50  0001 C CNN
F 3 "" H 2400 9000 50  0001 C CNN
	1    2400 9000
	0    1    1    0   
$EndComp
Text Label 2400 9200 2    50   ~ 0
3.3V
Text Label 3150 9200 0    50   ~ 0
3.3V
$Comp
L power:GND #PWR?
U 1 1 5F47759A
P 3150 9400
F 0 "#PWR?" H 3150 9150 50  0001 C CNN
F 1 "GND" V 3250 9350 50  0000 R CNN
F 2 "" H 3150 9400 50  0001 C CNN
F 3 "" H 3150 9400 50  0001 C CNN
	1    3150 9400
	0    -1   -1   0   
$EndComp
Text Label 3150 9000 0    50   ~ 0
A-D_IO_1
Text Label 2400 9400 2    50   ~ 0
A-D_IO_2
$Comp
L Custom_parts:BAT54SDW D?
U 1 1 5F48EF0C
P 2650 9800
F 0 "D?" H 2825 10147 60  0000 C CNN
F 1 "BAT54SDW" H 2825 10041 60  0000 C CNN
F 2 "Package_TO_SOT_SMD:SOT-363_SC-70-6" H 2850 10000 60  0001 L CNN
F 3 "https://www.diodes.com/assets/Datasheets/ds11005.pdf" H 2850 10100 60  0001 L CNN
F 4 "Diodes Incorporated" H 2650 9800 50  0001 C CNN "Manufacturer"
F 5 "BAT54SDW-7-F" H 2650 9800 50  0001 C CNN "Part #"
	1    2650 9800
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F48EF12
P 2450 9700
F 0 "#PWR?" H 2450 9450 50  0001 C CNN
F 1 "GND" V 2350 9650 50  0000 R CNN
F 2 "" H 2450 9700 50  0001 C CNN
F 3 "" H 2450 9700 50  0001 C CNN
	1    2450 9700
	0    1    1    0   
$EndComp
Text Label 2450 9900 2    50   ~ 0
3.3V
Text Label 3200 9900 0    50   ~ 0
3.3V
$Comp
L power:GND #PWR?
U 1 1 5F48EF1A
P 3200 10100
F 0 "#PWR?" H 3200 9850 50  0001 C CNN
F 1 "GND" V 3300 10050 50  0000 R CNN
F 2 "" H 3200 10100 50  0001 C CNN
F 3 "" H 3200 10100 50  0001 C CNN
	1    3200 10100
	0    -1   -1   0   
$EndComp
Text Label 3200 9700 0    50   ~ 0
A-D_IO_3
Text Label 2450 10100 2    50   ~ 0
A-D_IO_4
Text Notes 1800 8700 0    59   ~ 0
4-channel analog/digital IO with 0-3.3V clamp
$Comp
L Device:R_Pack08 RN?
U 1 1 5F4ABE10
P 2200 8150
F 0 "RN?" V 1583 8150 50  0000 C CNN
F 1 "10k" V 1674 8150 50  0000 C CNN
F 2 "Custom Footprints:Bourns_16-SOIC_8x_r-pack" V 2675 8150 50  0001 C CNN
F 3 "file://media.digikey.com/api/videos/videoplayer/smallplayer/N0706.pdf" H 2200 8150 50  0001 C CNN
F 4 "Bourns Inc." V 2200 8150 50  0001 C CNN "Manufacturer"
F 5 "4816P-T01-103LF" V 2200 8150 50  0001 C CNN "Part #"
	1    2200 8150
	0    1    1    0   
$EndComp
Wire Wire Line
	1550 7750 1550 8050
Wire Wire Line
	1550 8050 2000 8050
Wire Wire Line
	1550 8550 1550 8150
Wire Wire Line
	1550 8150 2000 8150
Wire Wire Line
	1550 8750 1600 8750
Wire Wire Line
	1600 8750 1600 8250
Wire Wire Line
	1600 8250 2000 8250
Wire Wire Line
	1550 9450 1650 9450
Wire Wire Line
	1650 9450 1650 8350
Wire Wire Line
	1650 8350 2000 8350
Wire Wire Line
	1550 9650 1700 9650
Wire Wire Line
	1700 9650 1700 8450
Wire Wire Line
	1700 8450 2000 8450
Wire Wire Line
	1550 7550 1600 7550
Wire Wire Line
	1600 7550 1600 7950
Wire Wire Line
	1600 7950 2000 7950
Wire Wire Line
	1550 6850 1650 6850
Wire Wire Line
	1650 6850 1650 7850
Wire Wire Line
	1650 7850 2000 7850
Wire Wire Line
	1550 6650 1700 6650
Wire Wire Line
	1700 6650 1700 7750
Wire Wire Line
	1700 7750 2000 7750
Text Label 2400 7750 0    50   ~ 0
external_analog_1
Text Label 2600 7300 2    50   ~ 0
external_analog_2
Text Label 3200 7300 0    50   ~ 0
external_analog_3
Text Label 3200 7100 0    50   ~ 0
external_analog_4
$Comp
L power:GND #PWR?
U 1 1 5F4EB97D
P 2600 7200
F 0 "#PWR?" H 2600 6950 50  0001 C CNN
F 1 "GND" V 2605 7072 50  0000 R CNN
F 2 "" H 2600 7200 50  0001 C CNN
F 3 "" H 2600 7200 50  0001 C CNN
	1    2600 7200
	0    1    1    0   
$EndComp
$Comp
L Custom_parts:D_Zener_x4_ACom_AKKKK D?
U 1 1 5F4F11C3
P 2900 7350
F 0 "D?" H 2900 6875 50  0000 C CNN
F 1 "D_Zener_x4_ACom_AKKKK" H 2900 6966 50  0000 C CNN
F 2 "Custom Footprints:SOT-753" H 2900 7100 50  0001 C CNN
F 3 "https://rohmfs.rohm.com/api/videos/videoplayer/smallplayer/ftz5.6e.pdf" H 2900 7100 50  0001 C CNN
F 4 "Rohm Semiconductor" H 2900 7350 50  0001 C CNN "Manufacturer"
F 5 "FTZ5.6ET148" H 2900 7350 50  0001 C CNN "Part #"
	1    2900 7350
	-1   0    0    1   
$EndComp
Text Notes 1800 6750 0    59   ~ 0
4-channel external analog input with 5.6V zener
Text Label 4900 7200 2    50   ~ 0
external_analog_1
Text Label 5500 7200 0    50   ~ 0
external_analog_2
Text Label 4900 7300 2    50   ~ 0
external_analog_2
Text Label 5500 7300 0    50   ~ 0
external_analog_3
Text Label 4900 7400 2    50   ~ 0
external_analog_3
Text Label 5500 7400 0    50   ~ 0
external_analog_4
$Comp
L Connector:Barrel_Jack_Switch J?
U 1 1 5F4FA1F9
P 1250 2100
F 0 "J?" H 1307 2417 50  0000 C CNN
F 1 "Barrel_Jack_Switch" H 1307 2326 50  0000 C CNN
F 2 "Custom Footprints:54-00165-DC_Jack" H 1300 2060 50  0001 C CNN
F 3 "https://www.tensility.com/api/videos/videoplayer/smallplayer/54-00165.pdf" H 1300 2060 50  0001 C CNN
F 4 "Tensility International Corp" H 1250 2100 50  0001 C CNN "Manufacturer"
F 5 "54-00165" H 1250 2100 50  0001 C CNN "Part #"
	1    1250 2100
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F5035F8
P 1550 2200
F 0 "#PWR?" H 1550 1950 50  0001 C CNN
F 1 "GND" H 1750 2100 50  0000 R CNN
F 2 "" H 1550 2200 50  0001 C CNN
F 3 "" H 1550 2200 50  0001 C CNN
	1    1550 2200
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Female J?
U 1 1 5F50593C
P 1200 1500
F 0 "J?" H 1050 1600 50  0000 L CNN
F 1 "Conn_01x02_Female" H 800 1350 50  0000 L CNN
F 2 "Connector_Wire:SolderWire-2.5sqmm_1x02_P7.2mm_D2.4mm_OD3.6mm" H 1200 1500 50  0001 C CNN
F 3 "~" H 1200 1500 50  0001 C CNN
	1    1200 1500
	-1   0    0    -1  
$EndComp
Wire Wire Line
	1550 2000 1800 2000
Wire Wire Line
	1550 2100 1550 2200
Connection ~ 1550 2200
Wire Wire Line
	1800 2000 1800 1500
Wire Wire Line
	1800 1500 1400 1500
Connection ~ 1800 2000
Wire Wire Line
	1800 2000 2200 2000
$Comp
L power:GND #PWR?
U 1 1 5F518EFF
P 1400 1600
F 0 "#PWR?" H 1400 1350 50  0001 C CNN
F 1 "GND" V 1405 1472 50  0000 R CNN
F 2 "" H 1400 1600 50  0001 C CNN
F 3 "" H 1400 1600 50  0001 C CNN
	1    1400 1600
	0    -1   -1   0   
$EndComp
$Comp
L Device:R_Pack08 RN?
U 1 1 5F51AB0E
P 2150 11500
F 0 "RN?" V 1533 11500 50  0000 C CNN
F 1 "10k" V 1624 11500 50  0000 C CNN
F 2 "Custom Footprints:Bourns_16-SOIC_8x_r-pack" V 2625 11500 50  0001 C CNN
F 3 "file://media.digikey.com/api/videos/videoplayer/smallplayer/N0706.pdf" H 2150 11500 50  0001 C CNN
F 4 "Bourns Inc." V 2150 11500 50  0001 C CNN "Manufacturer"
F 5 "4816P-T01-103LF" V 2150 11500 50  0001 C CNN "Part #"
	1    2150 11500
	0    1    1    0   
$EndComp
Text Label 2350 11100 0    50   ~ 0
Isense_1
Text Label 2350 11300 0    50   ~ 0
Isense_2
Text Label 2350 11500 0    50   ~ 0
Isense_3
Text Label 2350 11700 0    50   ~ 0
Isense_4
Wire Wire Line
	1950 11100 1950 11150
Wire Wire Line
	1950 11300 1950 11350
Wire Wire Line
	1950 11500 1950 11550
Wire Wire Line
	1950 11700 1950 11750
$Comp
L power:GND #PWR?
U 1 1 5F5306B0
P 2350 11200
F 0 "#PWR?" H 2350 10950 50  0001 C CNN
F 1 "GND" V 2350 11050 50  0000 R CNN
F 2 "" H 2350 11200 50  0001 C CNN
F 3 "" H 2350 11200 50  0001 C CNN
	1    2350 11200
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F531E2E
P 2350 11400
F 0 "#PWR?" H 2350 11150 50  0001 C CNN
F 1 "GND" V 2350 11250 50  0000 R CNN
F 2 "" H 2350 11400 50  0001 C CNN
F 3 "" H 2350 11400 50  0001 C CNN
	1    2350 11400
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F53235A
P 2350 11600
F 0 "#PWR?" H 2350 11350 50  0001 C CNN
F 1 "GND" V 2350 11450 50  0000 R CNN
F 2 "" H 2350 11600 50  0001 C CNN
F 3 "" H 2350 11600 50  0001 C CNN
	1    2350 11600
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F5328BE
P 2350 11800
F 0 "#PWR?" H 2350 11550 50  0001 C CNN
F 1 "GND" V 2350 11650 50  0000 R CNN
F 2 "" H 2350 11800 50  0001 C CNN
F 3 "" H 2350 11800 50  0001 C CNN
	1    2350 11800
	0    -1   -1   0   
$EndComp
$Comp
L Connector:Conn_Coaxial_x2 J?
U 1 1 5F532BE3
P 1450 11250
F 0 "J?" H 1400 11600 50  0000 L CNN
F 1 "Conn_Coaxial_x2" H 1150 11500 50  0000 L CNN
F 2 "Custom Footprints:031-6575_2x_BNC" H 1450 11150 50  0001 C CNN
F 3 "https://www.amphenolrf.com/library/download/link/link_id/436470/parent/031-6575/" H 1450 11150 50  0001 C CNN
F 4 "Amphenol RF" H 1450 11250 50  0001 C CNN "Manufacturer"
F 5 "031-6575" H 1450 11250 50  0001 C CNN "Part #"
	1    1450 11250
	-1   0    0    -1  
$EndComp
$Comp
L Connector:Conn_Coaxial_x2 J?
U 1 1 5F535BEA
P 1600 11650
F 0 "J?" H 1550 11900 50  0000 L CNN
F 1 "Conn_Coaxial_x2" H 950 11350 50  0000 L CNN
F 2 "Custom Footprints:031-6575_2x_BNC" H 1600 11550 50  0001 C CNN
F 3 "https://www.amphenolrf.com/library/download/link/link_id/436470/parent/031-6575/" H 1600 11550 50  0001 C CNN
F 4 "Amphenol RF" H 1600 11650 50  0001 C CNN "Manufacturer"
F 5 "031-6575" H 1600 11650 50  0001 C CNN "Part #"
	1    1600 11650
	-1   0    0    -1  
$EndComp
Wire Wire Line
	1650 11150 1950 11150
Connection ~ 1950 11150
Wire Wire Line
	1950 11150 1950 11200
Wire Wire Line
	1650 11350 1950 11350
Connection ~ 1950 11350
Wire Wire Line
	1950 11350 1950 11400
Wire Wire Line
	1800 11550 1950 11550
Connection ~ 1950 11550
Wire Wire Line
	1950 11550 1950 11600
Wire Wire Line
	1800 11750 1950 11750
Connection ~ 1950 11750
Wire Wire Line
	1950 11750 1950 11800
$Comp
L power:GND #PWR?
U 1 1 5F54BD71
P 1450 11550
F 0 "#PWR?" H 1450 11300 50  0001 C CNN
F 1 "GND" H 1550 11400 50  0000 R CNN
F 2 "" H 1450 11550 50  0001 C CNN
F 3 "" H 1450 11550 50  0001 C CNN
	1    1450 11550
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F54C869
P 1600 11950
F 0 "#PWR?" H 1600 11700 50  0001 C CNN
F 1 "GND" H 1700 11800 50  0000 R CNN
F 2 "" H 1600 11950 50  0001 C CNN
F 3 "" H 1600 11950 50  0001 C CNN
	1    1600 11950
	1    0    0    -1  
$EndComp
Wire Wire Line
	1800 2000 1800 2350
Wire Wire Line
	1800 2350 2050 2350
$EndSCHEMATC
