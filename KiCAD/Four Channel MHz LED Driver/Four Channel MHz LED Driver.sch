EESchema Schematic File Version 4
EELAYER 30 0
EELAYER END
$Descr A3 16535 11693
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
P 2700 7200
F 0 "U?" H 3350 6563 60  0000 C CNN
F 1 "DAC084S085" H 3350 6669 60  0000 C CNN
F 2 "Custom Footprints:DAC084S085CIMM" H 3400 7400 60  0001 C CNN
F 3 "http://www.ti.com/general/docs/suppproductinfo.tsp?distId=10&gotoUrl=http%3A%2F%2Fwww.ti.com%2Flit%2Fgpn%2Fdac084s085" H 3350 6669 60  0001 C CNN
F 4 "Texas Instruments" H 2700 7200 50  0001 C CNN "Manufacturer"
F 5 "DAC084S085CIMM/NOPB" H 2700 7200 50  0001 C CNN "Part #"
	1    2700 7200
	-1   0    0    1   
$EndComp
$Comp
L Custom_parts:MCP4351-502E_ST U?
U 1 1 5F3407BD
P 3900 6550
F 0 "U?" H 4500 6807 60  0000 C CNN
F 1 "MCP4211-502E_ST" H 4500 6701 60  0000 C CNN
F 2 "Custom Footprints:MCP4351-502E" H 4550 6800 60  0001 C CNN
F 3 "http://www.microchip.com/mymicrochip/filehandler.aspx?ddocname=en547555" H 4500 6701 60  0001 C CNN
F 4 "Microchip Technology" H 3900 6550 50  0001 C CNN "Manufacturer"
F 5 "MCP4351-502E/ST" H 3900 6550 50  0001 C CNN "Part #"
	1    3900 6550
	1    0    0    -1  
$EndComp
Wire Wire Line
	5100 7450 5100 7550
Wire Wire Line
	5100 6250 5100 6550
$Comp
L power:GND #PWR?
U 1 1 5F34784A
P 3900 6750
F 0 "#PWR?" H 3900 6500 50  0001 C CNN
F 1 "GND" V 3905 6622 50  0000 R CNN
F 2 "" H 3900 6750 50  0001 C CNN
F 3 "" H 3900 6750 50  0001 C CNN
	1    3900 6750
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F34877E
P 3900 7250
F 0 "#PWR?" H 3900 7000 50  0001 C CNN
F 1 "GND" V 3800 7250 50  0000 R CNN
F 2 "" H 3900 7250 50  0001 C CNN
F 3 "" H 3900 7250 50  0001 C CNN
	1    3900 7250
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F348E43
P 5100 7250
F 0 "#PWR?" H 5100 7000 50  0001 C CNN
F 1 "GND" V 5105 7122 50  0000 R CNN
F 2 "" H 5100 7250 50  0001 C CNN
F 3 "" H 5100 7250 50  0001 C CNN
	1    5100 7250
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F349663
P 5100 6750
F 0 "#PWR?" H 5100 6500 50  0001 C CNN
F 1 "GND" V 5105 6622 50  0000 R CNN
F 2 "" H 5100 6750 50  0001 C CNN
F 3 "" H 5100 6750 50  0001 C CNN
	1    5100 6750
	0    -1   -1   0   
$EndComp
Wire Wire Line
	3900 7150 3900 7250
Connection ~ 3900 7250
NoConn ~ 5100 6950
NoConn ~ 5100 7050
Text Label 5100 6850 0    50   ~ 0
5V-analog
$Comp
L power:GND #PWR?
U 1 1 5F34BA64
P 1400 6800
F 0 "#PWR?" H 1400 6550 50  0001 C CNN
F 1 "GND" V 1405 6672 50  0000 R CNN
F 2 "" H 1400 6800 50  0001 C CNN
F 3 "" H 1400 6800 50  0001 C CNN
	1    1400 6800
	0    1    1    0   
$EndComp
Text Label 1400 6900 2    50   ~ 0
5V-analog
$Comp
L power:GND #PWR?
U 1 1 5F3504EC
P 2700 7400
F 0 "#PWR?" H 2700 7150 50  0001 C CNN
F 1 "GND" H 2750 7250 50  0000 R CNN
F 2 "" H 2700 7400 50  0001 C CNN
F 3 "" H 2700 7400 50  0001 C CNN
	1    2700 7400
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F352F78
P 5200 7050
F 0 "#PWR?" H 5200 6800 50  0001 C CNN
F 1 "GND" H 5400 6950 50  0000 R CNN
F 2 "" H 5200 7050 50  0001 C CNN
F 3 "" H 5200 7050 50  0001 C CNN
	1    5200 7050
	1    0    0    -1  
$EndComp
Wire Wire Line
	5200 6850 5100 6850
$Comp
L Custom_parts:S18V20F12_12V_DC U?
U 1 1 5F357FC0
P 2350 1300
F 0 "U?" H 2325 1435 50  0000 C CNN
F 1 "S20V18F12_12V_DC" H 2325 1344 50  0000 C CNN
F 2 "Custom Footprints:S18V20F12_12V_DC" H 2350 1300 50  0001 C CNN
F 3 "https://www.pololu.com/product-info-merged/2577" H 2350 1300 50  0001 C CNN
F 4 "Pololu Corporation" H 2350 1300 50  0001 C CNN "Manufacturer"
F 5 "2577" H 2350 1300 50  0001 C CNN "Part #"
	1    2350 1300
	1    0    0    -1  
$EndComp
Wire Wire Line
	1850 1600 1850 1650
Wire Wire Line
	2800 1600 2800 1650
Text Label 1750 1650 2    50   ~ 0
Vin
Connection ~ 1850 1650
Wire Wire Line
	1850 1650 1850 1700
Wire Wire Line
	2800 1650 2900 1650
Connection ~ 2800 1650
Wire Wire Line
	2800 1650 2800 1700
Text Label 2900 1650 0    50   ~ 0
12V
$Comp
L Custom_parts:DPBW06F-05 U?
U 1 1 5F35D26E
P 4200 1300
F 0 "U?" H 4250 1435 50  0000 C CNN
F 1 "DPBW06F-05" H 4250 1344 50  0000 C CNN
F 2 "Custom Footprints:DPBW06F-05" H 4200 1300 50  0001 C CNN
F 3 "https://media.digikey.com/pdf/Data%20Sheets/Mean%20Well%20PDF's/SPBW06,DPBW06_Ds.pdf" H 4200 1300 50  0001 C CNN
F 4 "MEAN WELL USA Inc." H 4200 1300 50  0001 C CNN "Manufacturer"
F 5 "DPBW06F-05" H 4200 1300 50  0001 C CNN "Part #"
	1    4200 1300
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F35E7FC
P 2900 1500
F 0 "C?" H 2992 1546 50  0000 L CNN
F 1 "47uF" H 2992 1455 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 2900 1500 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 2900 1500 50  0001 C CNN
F 4 "Murata Electronics" H 2900 1500 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 2900 1500 50  0001 C CNN "Part #"
	1    2900 1500
	1    0    0    -1  
$EndComp
NoConn ~ 3750 1600
$Comp
L power:GND #PWR?
U 1 1 5F3602CD
P 3750 1400
F 0 "#PWR?" H 3750 1150 50  0001 C CNN
F 1 "GND" V 3755 1272 50  0000 R CNN
F 2 "" H 3750 1400 50  0001 C CNN
F 3 "" H 3750 1400 50  0001 C CNN
	1    3750 1400
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F36457B
P 5200 1500
F 0 "#PWR?" H 5200 1250 50  0001 C CNN
F 1 "GND" V 5205 1372 50  0000 R CNN
F 2 "" H 5200 1500 50  0001 C CNN
F 3 "" H 5200 1500 50  0001 C CNN
	1    5200 1500
	0    -1   -1   0   
$EndComp
$Comp
L pspice:INDUCTOR L?
U 1 1 5F3651AC
P 3300 1650
F 0 "L?" H 3300 1865 50  0000 C CNN
F 1 "10uH" H 3300 1774 50  0000 C CNN
F 2 "Inductor_SMD:L_1210_3225Metric" H 3300 1650 50  0001 C CNN
F 3 "https://product.tdk.com/info/en/catalog/datasheets/inductor_automotive_power_tfm322512alma_en.pdf?ref_disty=digikey" H 3300 1650 50  0001 C CNN
F 4 "TDK Corporation" H 3300 1650 50  0001 C CNN "Manufacturer"
F 5 "TFM322512ALMA100MTAA" H 3300 1650 50  0001 C CNN "Part #"
	1    3300 1650
	1    0    0    -1  
$EndComp
Wire Wire Line
	3050 1650 2900 1650
Connection ~ 2900 1650
Wire Wire Line
	3550 1650 3550 1500
Wire Wire Line
	3550 1500 3750 1500
$Comp
L Device:C_Small C?
U 1 1 5F36881C
P 4250 1850
F 0 "C?" V 4150 1850 50  0000 C CNN
F 1 "1500pF" V 4050 1850 50  0000 C CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 4250 1850 50  0001 C CNN
F 3 "https://katalog.we-online.de/pbs/datasheet/885342206004.pdf" H 4250 1850 50  0001 C CNN
F 4 "Würth Elektronik" H 4250 1850 50  0001 C CNN "Manufacturer"
F 5 "885342206004" H 4250 1850 50  0001 C CNN "Part #"
	1    4250 1850
	0    1    -1   0   
$EndComp
Wire Wire Line
	4150 1850 3550 1850
Wire Wire Line
	3550 1850 3550 1650
Connection ~ 3550 1650
Wire Wire Line
	4350 1850 4750 1850
$Comp
L Device:C_Small C?
U 1 1 5F36C3BB
P 4250 1050
F 0 "C?" V 4350 1050 50  0000 C CNN
F 1 "1500pF" V 4300 1250 50  0000 C CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 4250 1050 50  0001 C CNN
F 3 "https://katalog.we-online.de/pbs/datasheet/885342206004.pdf" H 4250 1050 50  0001 C CNN
F 4 "Würth Elektronik" H 4250 1050 50  0001 C CNN "Manufacturer"
F 5 "885342206004" H 4250 1050 50  0001 C CNN "Part #"
	1    4250 1050
	0    1    -1   0   
$EndComp
Wire Wire Line
	3750 1400 3750 1050
Wire Wire Line
	3750 1050 4150 1050
Connection ~ 3750 1400
Wire Wire Line
	4750 1400 4750 1300
Wire Wire Line
	4750 1050 4350 1050
$Comp
L Device:C_Small C?
U 1 1 5F370F72
P 4950 1400
F 0 "C?" H 5042 1446 50  0000 L CNN
F 1 "47uF" H 5042 1355 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 4950 1400 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 4950 1400 50  0001 C CNN
F 4 "Murata Electronics" H 4950 1400 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 4950 1400 50  0001 C CNN "Part #"
	1    4950 1400
	1    0    0    -1  
$EndComp
Connection ~ 4750 1300
Wire Wire Line
	4750 1300 4750 1050
Text Label 4800 1700 0    50   ~ 0
5V
Text Label 4750 1300 0    50   ~ 0
-5V
Connection ~ 4950 1500
Wire Wire Line
	4750 1500 4950 1500
Wire Wire Line
	4750 1300 4950 1300
Wire Wire Line
	2900 1600 2900 1650
Wire Wire Line
	2800 1400 2800 1500
Wire Wire Line
	2800 1400 2900 1400
Connection ~ 2800 1400
$Comp
L Device:R_Small R?
U 1 1 5F379C80
P 2900 1750
F 0 "R?" H 2959 1796 50  0000 L CNN
F 1 "100" H 2959 1705 50  0000 L CNN
F 2 "Resistor_SMD:R_0805_2012Metric" H 2900 1750 50  0001 C CNN
F 3 "https://www.seielect.com/catalog/sei-rncp.pdf" H 2900 1750 50  0001 C CNN
F 4 "Stackpole Electronics Inc" H 2900 1750 50  0001 C CNN "Manufacturer"
F 5 "RNCP0805FTD100R" H 2900 1750 50  0001 C CNN "Part #"
	1    2900 1750
	1    0    0    -1  
$EndComp
$Comp
L Custom_parts:ADP7118AUJZ-5.0-R7 U?
U 1 1 5F380EEB
P 3200 2350
F 0 "U?" H 3625 2557 60  0000 C CNN
F 1 "ADP7120AUJZ-5.0-R7" H 3625 2451 60  0000 C CNN
F 2 "Custom Footprints:ADP7118AUJZ-R7" H 4300 2590 60  0001 C CNN
F 3 "https://www.analog.com/media/en/technical-documentation/data-sheets/ADP7118.pdf" H 3625 2451 60  0001 C CNN
F 4 "Analog Devices Inc." H 3200 2350 50  0001 C CNN "Manufacturer"
F 5 "ADP7118AUJZ-5.0-R7" H 3200 2350 50  0001 C CNN "Part #"
	1    3200 2350
	1    0    0    -1  
$EndComp
Connection ~ 2900 2400
Wire Wire Line
	4350 2600 4350 2400
Text Label 4350 2500 0    50   ~ 0
5V-analog
Wire Wire Line
	2900 1850 2900 2000
$Comp
L Device:C_Small C?
U 1 1 5F38C930
P 2800 2400
F 0 "C?" V 3000 2350 50  0000 L CNN
F 1 "2.2uF" V 2900 2300 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 2800 2400 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 2800 2400 50  0001 C CNN
F 4 "Taiyo Yuden" H 2800 2400 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 2800 2400 50  0001 C CNN "Part #"
	1    2800 2400
	0    -1   -1   0   
$EndComp
Connection ~ 2900 2000
$Comp
L power:GND #PWR?
U 1 1 5F39026E
P 2700 2000
F 0 "#PWR?" H 2700 1750 50  0001 C CNN
F 1 "GND" V 2700 1850 50  0000 R CNN
F 2 "" H 2700 2000 50  0001 C CNN
F 3 "" H 2700 2000 50  0001 C CNN
	1    2700 2000
	0    1    1    0   
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F390986
P 5200 6950
F 0 "C?" H 5292 6996 50  0000 L CNN
F 1 "2.2uF" H 5292 6905 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 5200 6950 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 5200 6950 50  0001 C CNN
F 4 "Taiyo Yuden" H 5200 6950 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 5200 6950 50  0001 C CNN "Part #"
	1    5200 6950
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F3918D7
P 2700 7300
F 0 "C?" H 2500 7300 50  0000 L CNN
F 1 "2.2uF" H 2450 7200 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 2700 7300 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 2700 7300 50  0001 C CNN
F 4 "Taiyo Yuden" H 2700 7300 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 2700 7300 50  0001 C CNN "Part #"
	1    2700 7300
	1    0    0    -1  
$EndComp
Text Label 2700 7200 0    50   ~ 0
5V-analog
Wire Wire Line
	2700 7000 3200 7000
Wire Wire Line
	2700 7100 3150 7100
Text Notes 2450 2750 0    59   ~ 0
LDO - 12V to clean 5V for analog circuits
Text Notes 3450 900  0    59   ~ 0
ISOLATED - 12V to split +/- 5V
Text Notes 1550 650  0    59   ~ 0
SEPIC - Vin (3V - 30V) to 12V DC
Wire Wire Line
	2900 2000 2900 2100
$Comp
L Jumper:SolderJumper_2_Open JP?
U 1 1 5F3A0F69
P 1850 2000
F 0 "JP?" H 1850 2100 50  0000 C CNN
F 1 "SolderJumper_2_Open" H 1850 1900 50  0000 C CNN
F 2 "Jumper:SolderJumper-2_P1.3mm_Open_Pad1.0x1.5mm" H 1850 2000 50  0001 C CNN
F 3 "~" H 1850 2000 50  0001 C CNN
	1    1850 2000
	1    0    0    -1  
$EndComp
Text Label 2000 2000 0    50   ~ 0
12V
$Comp
L Device:C_Small C?
U 1 1 5F3A29AA
P 2800 2000
F 0 "C?" V 2571 2000 50  0000 C CNN
F 1 "47uF" V 2662 2000 50  0000 C CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 2800 2000 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 2800 2000 50  0001 C CNN
F 4 "Murata Electronics" H 2800 2000 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 2800 2000 50  0001 C CNN "Part #"
	1    2800 2000
	0    1    1    0   
$EndComp
Connection ~ 2900 2100
Wire Wire Line
	2900 2100 2900 2400
$Comp
L power:GND #PWR?
U 1 1 5F3A55D8
P 2700 2400
F 0 "#PWR?" H 2700 2150 50  0001 C CNN
F 1 "GND" V 2700 2250 50  0000 R CNN
F 2 "" H 2700 2400 50  0001 C CNN
F 3 "" H 2700 2400 50  0001 C CNN
	1    2700 2400
	0    1    1    0   
$EndComp
Wire Wire Line
	2700 2400 2700 2500
Wire Wire Line
	2700 2500 2900 2500
Connection ~ 2700 2400
Wire Wire Line
	2350 2600 2350 2100
Wire Wire Line
	2350 2600 2900 2600
Wire Wire Line
	2350 2100 2900 2100
Text Notes 1400 2250 0    59   ~ 0
12V Supply Bypass
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3BA201
P 6750 7300
F 0 "U?" H 6925 7465 50  0000 C CNN
F 1 "TMUX1184DGSR" H 6925 7374 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 6750 8300 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 6750 7300 50  0001 C CNN
F 4 "Texas Instruments" H 6750 7300 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 6750 7300 50  0001 C CNN "Part #"
	1    6750 7300
	1    0    0    -1  
$EndComp
Wire Wire Line
	6550 7750 6550 7850
Wire Wire Line
	6550 7850 7300 7850
Wire Wire Line
	7300 7850 7300 7750
Wire Wire Line
	7300 7450 7350 7450
Wire Wire Line
	7350 7450 7350 7650
Wire Wire Line
	7350 7650 7300 7650
$Comp
L power:GND #PWR?
U 1 1 5F3BE358
P 6550 7550
F 0 "#PWR?" H 6550 7300 50  0001 C CNN
F 1 "GND" V 6555 7422 50  0000 R CNN
F 2 "" H 6550 7550 50  0001 C CNN
F 3 "" H 6550 7550 50  0001 C CNN
	1    6550 7550
	0    1    1    0   
$EndComp
Text Label 5100 7350 0    50   ~ 0
internal_analog_1
Text Label 5100 6650 0    50   ~ 0
internal_analog_3
Text Label 3900 6650 2    50   ~ 0
internal_analog_4
Text Label 3900 7350 2    50   ~ 0
internal_analog_2
Wire Wire Line
	3200 7450 3200 7000
Wire Wire Line
	3200 7450 3900 7450
Wire Wire Line
	3150 7100 3150 7550
Wire Wire Line
	3150 7550 5100 7550
Text Label 6550 7450 2    50   ~ 0
internal_analog_1
Text Label 6550 7650 2    50   ~ 0
external_analog_1
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3DA347
P 6750 8150
F 0 "U?" H 6925 8315 50  0000 C CNN
F 1 "TMUX1184DGSR" H 6925 8224 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 6750 9150 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 6750 8150 50  0001 C CNN
F 4 "Texas Instruments" H 6750 8150 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 6750 8150 50  0001 C CNN "Part #"
	1    6750 8150
	1    0    0    -1  
$EndComp
Wire Wire Line
	6550 8600 6550 8700
Wire Wire Line
	6550 8700 7300 8700
Wire Wire Line
	7300 8700 7300 8600
Wire Wire Line
	7300 8300 7350 8300
Wire Wire Line
	7350 8300 7350 8500
Wire Wire Line
	7350 8500 7300 8500
$Comp
L power:GND #PWR?
U 1 1 5F3DA355
P 6550 8400
F 0 "#PWR?" H 6550 8150 50  0001 C CNN
F 1 "GND" V 6555 8272 50  0000 R CNN
F 2 "" H 6550 8400 50  0001 C CNN
F 3 "" H 6550 8400 50  0001 C CNN
	1    6550 8400
	0    1    1    0   
$EndComp
Text Label 6550 8300 2    50   ~ 0
internal_analog_2
Text Label 6550 8500 2    50   ~ 0
external_analog_2
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3DC268
P 6800 9000
F 0 "U?" H 6975 9165 50  0000 C CNN
F 1 "TMUX1184DGSR" H 6975 9074 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 6800 10000 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 6800 9000 50  0001 C CNN
F 4 "Texas Instruments" H 6800 9000 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 6800 9000 50  0001 C CNN "Part #"
	1    6800 9000
	1    0    0    -1  
$EndComp
Wire Wire Line
	6600 9450 6600 9550
Wire Wire Line
	6600 9550 7350 9550
Wire Wire Line
	7350 9550 7350 9450
Wire Wire Line
	7350 9150 7400 9150
Wire Wire Line
	7400 9150 7400 9350
Wire Wire Line
	7400 9350 7350 9350
$Comp
L power:GND #PWR?
U 1 1 5F3DC276
P 6600 9250
F 0 "#PWR?" H 6600 9000 50  0001 C CNN
F 1 "GND" V 6605 9122 50  0000 R CNN
F 2 "" H 6600 9250 50  0001 C CNN
F 3 "" H 6600 9250 50  0001 C CNN
	1    6600 9250
	0    1    1    0   
$EndComp
Text Label 6600 9150 2    50   ~ 0
internal_analog_3
Text Label 6600 9350 2    50   ~ 0
external_analog_3
$Comp
L Custom_parts:TMUX1204DGSR U?
U 1 1 5F3DEF78
P 6800 9800
F 0 "U?" H 6975 9965 50  0000 C CNN
F 1 "TMUX1184DGSR" H 6975 9874 50  0000 C CNN
F 2 "Custom Footprints:TMUX1204DGSR" H 6800 10800 50  0001 L BNN
F 3 "https://www.ti.com/api/videos/videoplayer/smallplayer/suppproductinfo.tsp" H 6800 9800 50  0001 C CNN
F 4 "Texas Instruments" H 6800 9800 50  0001 C CNN "Manufacturer"
F 5 "TMUX1204DGSR" H 6800 9800 50  0001 C CNN "Part #"
	1    6800 9800
	1    0    0    -1  
$EndComp
Wire Wire Line
	6600 10250 6600 10350
Wire Wire Line
	6600 10350 7350 10350
Wire Wire Line
	7350 10350 7350 10250
Wire Wire Line
	7350 9950 7400 9950
Wire Wire Line
	7400 9950 7400 10150
Wire Wire Line
	7400 10150 7350 10150
$Comp
L power:GND #PWR?
U 1 1 5F3DEF86
P 6600 10050
F 0 "#PWR?" H 6600 9800 50  0001 C CNN
F 1 "GND" V 6605 9922 50  0000 R CNN
F 2 "" H 6600 10050 50  0001 C CNN
F 3 "" H 6600 10050 50  0001 C CNN
	1    6600 10050
	0    1    1    0   
$EndComp
Text Label 6600 9950 2    50   ~ 0
internal_analog_4
Text Label 6600 10150 2    50   ~ 0
external_analog_4
Text Notes 1950 6100 0    59   ~ 0
Internal 4-channel ADC AWG with digipot current limit
Text Label 4250 8150 2    50   ~ 0
external_analog_1
Text Label 2650 8400 0    50   ~ 0
external_analog_2
Text Label 2650 8500 0    50   ~ 0
external_analog_3
Text Label 2650 8600 0    50   ~ 0
external_analog_4
$Comp
L Switch:SW_DIP_x03 SW?
U 1 1 5F411870
P 4700 9200
F 0 "SW?" H 4700 9667 50  0000 C CNN
F 1 "SW_DIP_x03" H 4700 9576 50  0000 C CNN
F 2 "Custom Footprints:SW_DS04-254-2-03BK-SMT" H 4700 9200 50  0001 C CNN
F 3 "https://www.cuidevices.com/api/videos/videoplayer/smallplayer/ds04-254-smt.pdf" H 4700 9200 50  0001 C CNN
F 4 "CUI Devices" H 4700 9200 50  0001 C CNN "Manufacturer"
F 5 "DS04-254-2-03BK-SMT" H 4700 9200 50  0001 C CNN "Part #"
	1    4700 9200
	1    0    0    -1  
$EndComp
Connection ~ 4750 1700
Wire Wire Line
	4750 1700 4750 1600
Wire Wire Line
	4750 1850 4750 1700
Wire Wire Line
	4950 1500 5200 1500
Text Label 5650 1800 2    50   ~ 0
-0.25V_analog
$Comp
L power:GND #PWR?
U 1 1 5F3D1FD7
P 5650 2050
F 0 "#PWR?" H 5650 1800 50  0001 C CNN
F 1 "GND" V 5750 2000 50  0000 R CNN
F 2 "" H 5650 2050 50  0001 C CNN
F 3 "" H 5650 2050 50  0001 C CNN
	1    5650 2050
	1    0    0    -1  
$EndComp
$Comp
L Device:R_Small R?
U 1 1 5F3D1763
P 5650 1950
F 0 "R?" H 5709 1996 50  0000 L CNN
F 1 "249" H 5709 1905 50  0000 L CNN
F 2 "Resistor_SMD:R_0603_1608Metric" H 5650 1950 50  0001 C CNN
F 3 "https://www.seielect.com/api/videos/videoplayer/smallplayer/sei-rncp.pdf" H 5650 1950 50  0001 C CNN
F 4 "Stackpole Electronics Inc" H 5650 1950 50  0001 C CNN "Manufacturer"
F 5 "RNCP0603FTD249R" H 5650 1950 50  0001 C CNN "Part #"
	1    5650 1950
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F3D013F
P 5850 1600
F 0 "#PWR?" H 5850 1350 50  0001 C CNN
F 1 "GND" V 5950 1550 50  0000 R CNN
F 2 "" H 5850 1600 50  0001 C CNN
F 3 "" H 5850 1600 50  0001 C CNN
	1    5850 1600
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5650 1600 5650 1850
Connection ~ 5650 1600
Wire Wire Line
	5650 1500 5650 1600
$Comp
L Device:C_Small C?
U 1 1 5F3CD242
P 5750 1600
F 0 "C?" V 5850 1550 50  0000 L CNN
F 1 "47uF" V 5950 1550 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 5750 1600 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 5750 1600 50  0001 C CNN
F 4 "Murata Electronics" H 5750 1600 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 5750 1600 50  0001 C CNN "Part #"
	1    5750 1600
	0    1    1    0   
$EndComp
$Comp
L Device:R_Small R?
U 1 1 5F3CC1B4
P 5650 1400
F 0 "R?" H 5709 1446 50  0000 L CNN
F 1 "4.7k" H 5709 1355 50  0000 L CNN
F 2 "Resistor_SMD:R_0805_2012Metric" H 5650 1400 50  0001 C CNN
F 3 "https://www.seielect.com/catalog/sei-rncp.pdf" H 5650 1400 50  0001 C CNN
F 4 "Yageo" H 5650 1400 50  0001 C CNN "Manufacturer"
F 5 "RT0805FRE074K7L" H 5650 1400 50  0001 C CNN "Part #"
	1    5650 1400
	1    0    0    -1  
$EndComp
Wire Wire Line
	4750 1700 4950 1700
$Comp
L Device:C_Small C?
U 1 1 5F37216A
P 4950 1600
F 0 "C?" H 5042 1646 50  0000 L CNN
F 1 "47uF" H 5042 1555 50  0000 L CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 4950 1600 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 4950 1600 50  0001 C CNN
F 4 "Murata Electronics" H 4950 1600 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 4950 1600 50  0001 C CNN "Part #"
	1    4950 1600
	1    0    0    -1  
$EndComp
Wire Wire Line
	4950 1300 5650 1300
Connection ~ 4950 1300
Text Label 7350 7450 0    50   ~ 0
-0.25V_analog
Text Label 7350 8300 0    50   ~ 0
-0.25V_analog
Text Label 7400 9150 0    50   ~ 0
-0.25V_analog
Text Label 7400 9950 0    50   ~ 0
-0.25V_analog
Text Label 7300 7850 0    50   ~ 0
5V-analog
Text Label 7300 8700 0    50   ~ 0
5V-analog
Text Label 7350 9550 0    50   ~ 0
5V-analog
Text Label 7350 10350 0    50   ~ 0
5V-analog
$Comp
L Device:R_Pack08 RN?
U 1 1 5F4ABE10
P 2450 8700
F 0 "RN?" V 1833 8700 50  0000 C CNN
F 1 "4.7k" V 1924 8700 50  0000 C CNN
F 2 "Custom Footprints:Bourns_16-SOIC_8x_r-pack" V 2925 8700 50  0001 C CNN
F 3 "file://media.digikey.com/api/videos/videoplayer/smallplayer/N0706.pdf" H 2450 8700 50  0001 C CNN
F 4 "Bourns Inc." V 2450 8700 50  0001 C CNN "Manufacturer"
F 5 "4816P-T01-472LF" V 2450 8700 50  0001 C CNN "Part #"
	1    2450 8700
	0    1    1    0   
$EndComp
Text Label 2650 8300 0    50   ~ 0
external_analog_1
Text Label 4250 8350 2    50   ~ 0
external_analog_2
Text Label 4850 8350 0    50   ~ 0
external_analog_3
Text Label 4850 8150 0    50   ~ 0
external_analog_4
$Comp
L power:GND #PWR?
U 1 1 5F4EB97D
P 4250 8250
F 0 "#PWR?" H 4250 8000 50  0001 C CNN
F 1 "GND" V 4255 8122 50  0000 R CNN
F 2 "" H 4250 8250 50  0001 C CNN
F 3 "" H 4250 8250 50  0001 C CNN
	1    4250 8250
	0    1    1    0   
$EndComp
$Comp
L Custom_parts:D_Zener_x4_ACom_AKKKK D?
U 1 1 5F4F11C3
P 4550 8400
F 0 "D?" H 4550 7925 50  0000 C CNN
F 1 "D_Zener_x4_ACom_AKKKK" H 4550 8016 50  0000 C CNN
F 2 "Custom Footprints:SOT-753" H 4550 8150 50  0001 C CNN
F 3 "https://rohmfs.rohm.com/api/videos/videoplayer/smallplayer/ftz5.6e.pdf" H 4550 8150 50  0001 C CNN
F 4 "Rohm Semiconductor" H 4550 8400 50  0001 C CNN "Manufacturer"
F 5 "FTZ5.6ET148" H 4550 8400 50  0001 C CNN "Part #"
	1    4550 8400
	-1   0    0    1   
$EndComp
Text Notes 3450 7800 0    59   ~ 0
4-channel external analog input with 5.6V zener
Text Label 4400 9000 2    50   ~ 0
external_analog_1
Text Label 5000 9000 0    50   ~ 0
external_analog_2
Text Label 4400 9100 2    50   ~ 0
external_analog_2
Text Label 5000 9100 0    50   ~ 0
external_analog_3
Text Label 4400 9200 2    50   ~ 0
external_analog_3
Text Label 5000 9200 0    50   ~ 0
external_analog_4
$Comp
L Connector:Barrel_Jack_Switch J?
U 1 1 5F4FA1F9
P 900 1750
F 0 "J?" H 957 2067 50  0000 C CNN
F 1 "Barrel_Jack_Switch" H 957 1976 50  0000 C CNN
F 2 "Custom Footprints:54-00165-DC_Jack" H 950 1710 50  0001 C CNN
F 3 "https://www.tensility.com/api/videos/videoplayer/smallplayer/54-00165.pdf" H 950 1710 50  0001 C CNN
F 4 "Tensility International Corp" H 900 1750 50  0001 C CNN "Manufacturer"
F 5 "54-00165" H 900 1750 50  0001 C CNN "Part #"
	1    900  1750
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F5035F8
P 1200 1850
F 0 "#PWR?" H 1200 1600 50  0001 C CNN
F 1 "GND" H 1400 1750 50  0000 R CNN
F 2 "" H 1200 1850 50  0001 C CNN
F 3 "" H 1200 1850 50  0001 C CNN
	1    1200 1850
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Female J?
U 1 1 5F50593C
P 850 1150
F 0 "J?" H 700 1250 50  0000 L CNN
F 1 "Conn_01x02_Female" H 450 1000 50  0000 L CNN
F 2 "Connector_Wire:SolderWire-2.5sqmm_1x02_P7.2mm_D2.4mm_OD3.6mm" H 850 1150 50  0001 C CNN
F 3 "~" H 850 1150 50  0001 C CNN
	1    850  1150
	-1   0    0    -1  
$EndComp
Wire Wire Line
	1200 1750 1200 1850
Connection ~ 1200 1850
$Comp
L power:GND #PWR?
U 1 1 5F518EFF
P 1050 1250
F 0 "#PWR?" H 1050 1000 50  0001 C CNN
F 1 "GND" V 1055 1122 50  0000 R CNN
F 2 "" H 1050 1250 50  0001 C CNN
F 3 "" H 1050 1250 50  0001 C CNN
	1    1050 1250
	0    -1   -1   0   
$EndComp
Text Label 2050 10350 2    50   ~ 0
Isense_1
Text Label 2050 10550 2    50   ~ 0
Isense_2
Text Label 2050 10750 2    50   ~ 0
Isense_3
Text Label 2050 10950 2    50   ~ 0
Isense_4
$Comp
L power:GND #PWR?
U 1 1 5F5306B0
P 2250 11050
F 0 "#PWR?" H 2250 10800 50  0001 C CNN
F 1 "GND" V 2250 10900 50  0000 R CNN
F 2 "" H 2250 11050 50  0001 C CNN
F 3 "" H 2250 11050 50  0001 C CNN
	1    2250 11050
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F531E2E
P 2250 10850
F 0 "#PWR?" H 2250 10600 50  0001 C CNN
F 1 "GND" V 2250 10700 50  0000 R CNN
F 2 "" H 2250 10850 50  0001 C CNN
F 3 "" H 2250 10850 50  0001 C CNN
	1    2250 10850
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F53235A
P 2250 10650
F 0 "#PWR?" H 2250 10400 50  0001 C CNN
F 1 "GND" V 2250 10500 50  0000 R CNN
F 2 "" H 2250 10650 50  0001 C CNN
F 3 "" H 2250 10650 50  0001 C CNN
	1    2250 10650
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F5328BE
P 2250 10450
F 0 "#PWR?" H 2250 10200 50  0001 C CNN
F 1 "GND" V 2250 10300 50  0000 R CNN
F 2 "" H 2250 10450 50  0001 C CNN
F 3 "" H 2250 10450 50  0001 C CNN
	1    2250 10450
	0    1    1    0   
$EndComp
Text Label 14700 900  0    50   ~ 0
5V
Text Label 12400 2300 2    50   ~ 0
3.3V
Text Label 14700 1100 0    50   ~ 0
3.3V
$Comp
L power:GND #PWR?
U 1 1 5F3EF0B8
P 12400 900
F 0 "#PWR?" H 12400 650 50  0001 C CNN
F 1 "GND" V 12400 750 50  0000 R CNN
F 2 "" H 12400 900 50  0001 C CNN
F 3 "" H 12400 900 50  0001 C CNN
	1    12400 900 
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F3DEA30
P 14700 1000
F 0 "#PWR?" H 14700 750 50  0001 C CNN
F 1 "GND" V 14700 850 50  0000 R CNN
F 2 "" H 14700 1000 50  0001 C CNN
F 3 "" H 14700 1000 50  0001 C CNN
	1    14700 1000
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F3F47A9
P 14700 2300
F 0 "#PWR?" H 14700 2050 50  0001 C CNN
F 1 "GND" V 14700 2150 50  0000 R CNN
F 2 "" H 14700 2300 50  0001 C CNN
F 3 "" H 14700 2300 50  0001 C CNN
	1    14700 2300
	0    -1   -1   0   
$EndComp
$Comp
L Custom_parts:Teensy3.6 U?
U 1 1 5F413A64
P 13550 3050
F 0 "U?" H 13550 5487 60  0000 C CNN
F 1 "Teensy3.6" H 13550 5381 60  0000 C CNN
F 2 "" H 13550 3100 60  0001 C CNN
F 3 "https://media.digikey.com/pdf/Data%20Sheets/Sparkfun%20PDFs/DEV-14058_Web.pdf" H 13550 5381 60  0001 C CNN
F 4 "SparkFun Electronics" H 13550 3050 50  0001 C CNN "Manufacturer"
F 5 "DEV-14058" H 13550 3050 50  0001 C CNN "Part #"
	1    13550 3050
	1    0    0    -1  
$EndComp
NoConn ~ 12400 3300
NoConn ~ 12400 3400
NoConn ~ 12400 3500
NoConn ~ 12400 3600
NoConn ~ 12400 3700
NoConn ~ 12400 3800
NoConn ~ 12400 3900
NoConn ~ 12400 4000
NoConn ~ 12400 4100
NoConn ~ 12400 4250
NoConn ~ 12400 4350
NoConn ~ 12400 4450
NoConn ~ 12400 4550
NoConn ~ 12400 4650
NoConn ~ 12400 4750
NoConn ~ 12400 4850
NoConn ~ 12400 4950
NoConn ~ 12400 5050
NoConn ~ 12400 5150
NoConn ~ 14700 3350
NoConn ~ 14700 3450
NoConn ~ 14700 3550
NoConn ~ 14700 3650
NoConn ~ 14700 3750
NoConn ~ 14700 3850
NoConn ~ 14700 3950
NoConn ~ 14700 4050
NoConn ~ 14700 4150
NoConn ~ 14700 4250
NoConn ~ 14700 4350
NoConn ~ 14700 4450
NoConn ~ 14700 4550
NoConn ~ 14700 4650
NoConn ~ 14700 4750
NoConn ~ 14700 4850
NoConn ~ 14700 4950
NoConn ~ 14700 5050
NoConn ~ 14700 5150
Wire Wire Line
	1850 1450 1850 1500
Wire Wire Line
	1850 1400 1850 1450
Connection ~ 1850 1450
$Comp
L power:GND #PWR?
U 1 1 5F35A16A
P 1850 1450
F 0 "#PWR?" H 1850 1200 50  0001 C CNN
F 1 "GND" V 1855 1322 50  0000 R CNN
F 2 "" H 1850 1450 50  0001 C CNN
F 3 "" H 1850 1450 50  0001 C CNN
	1    1850 1450
	0    1    1    0   
$EndComp
NoConn ~ 1850 1250
Wire Wire Line
	1850 1150 1850 1050
Wire Wire Line
	1850 950  1850 900 
Text Label 1850 1150 0    50   ~ 0
Vin
$Comp
L power:GND #PWR?
U 1 1 5F4F5E01
P 2800 900
F 0 "#PWR?" H 2800 650 50  0001 C CNN
F 1 "GND" V 2900 850 50  0000 R CNN
F 2 "" H 2800 900 50  0001 C CNN
F 3 "" H 2800 900 50  0001 C CNN
	1    2800 900 
	0    1    1    0   
$EndComp
Connection ~ 1850 900 
Wire Wire Line
	1850 900  1850 850 
$Comp
L Custom_parts:Conn_01x05_Male J?
U 1 1 5F515681
P 1650 1050
F 0 "J?" H 1750 1400 50  0000 C CNN
F 1 "Conn_01x05_Male" H 1750 1300 50  0000 C CNN
F 2 "Custom Footprints:Ref_only" H 1650 1050 50  0001 C CNN
F 3 "http://suddendocs.samtec.com/catalog_english/tsm.pdf" H 1650 1050 50  0001 C CNN
F 4 "Samtec Inc." H 1650 1050 50  0001 C CNN "Manufacturer"
F 5 "TSM-105-01-T-SV" H 1650 1050 50  0001 C CNN "Part #"
	1    1650 1050
	1    0    0    -1  
$EndComp
$Comp
L Custom_parts:Conn_01x04_Male J?
U 1 1 5F551B78
P 3000 1050
F 0 "J?" H 3100 700 50  0000 R CNN
F 1 "Conn_01x04_Male" H 3400 800 50  0000 R CNN
F 2 "Custom Footprints:Ref_only" H 3000 1050 50  0001 C CNN
F 3 "http://suddendocs.samtec.com/catalog_english/tsm.pdf" H 3000 1050 50  0001 C CNN
F 4 "Samtec Inc." H 3000 1050 50  0001 C CNN "Manufacturer"
F 5 "TSM-104-01-T-SV" H 3000 1050 50  0001 C CNN "Part #"
	1    3000 1050
	-1   0    0    1   
$EndComp
Wire Wire Line
	2800 950  2800 900 
Wire Wire Line
	2800 1050 2800 1150
Text Label 2800 1150 2    50   ~ 0
12V
Connection ~ 2800 900 
Wire Wire Line
	2800 900  2800 850 
$Comp
L power:GND #PWR?
U 1 1 5F56462A
P 1850 900
F 0 "#PWR?" H 1850 650 50  0001 C CNN
F 1 "GND" V 1750 850 50  0000 R CNN
F 2 "" H 1850 900 50  0001 C CNN
F 3 "" H 1850 900 50  0001 C CNN
	1    1850 900 
	0    -1   -1   0   
$EndComp
Text Notes 2100 1100 0    50   ~ 0
Header pins\nto connect\n12V SEPIC \nDC-DC \nconverter
NoConn ~ 1850 1800
Text Label 7400 7550 0    50   ~ 0
OA1_input
Text Label 7400 8400 0    50   ~ 0
OA2_input
Text Label 7450 9250 0    50   ~ 0
OA3_input
Text Label 7450 10050 0    50   ~ 0
OA4_input
Wire Wire Line
	7400 7550 7300 7550
Wire Wire Line
	7400 8400 7300 8400
Wire Wire Line
	7450 9250 7350 9250
Wire Wire Line
	7450 10050 7350 10050
Text Notes 1650 10050 0    59   ~ 0
Current sense voltage ouput
Text Notes 6500 7050 0    59   ~ 0
Op-amp input mux
Text Label 12400 2700 2    50   ~ 0
SCK0
Text Label 12400 2800 2    50   ~ 0
MOSI0
Text Label 3900 6950 2    50   ~ 0
SCK0
Text Label 3900 7050 2    50   ~ 0
MOSI0
Text Label 1400 7200 2    50   ~ 0
SCK0
Text Label 1400 7000 2    50   ~ 0
MOSI0
Text Label 12400 2400 2    50   ~ 0
24
Text Label 12400 2500 2    50   ~ 0
25
Text Label 1400 7100 2    50   ~ 0
25
Text Label 3900 6850 2    50   ~ 0
24
Text Notes 11850 2400 0    50   ~ 0
Digipot CS
Text Notes 12000 2500 0    50   ~ 0
DAC CS
Text Notes 11800 3100 0    50   ~ 0
Isense_1
Text Notes 11800 3200 0    50   ~ 0
Isense_2
Text Notes 15000 3200 0    50   ~ 0
Isense_3
Text Notes 15000 3100 0    50   ~ 0
Isense_4
Text Label 12400 1200 2    50   ~ 0
2
Text Label 12400 1300 2    50   ~ 0
3
Text Label 12400 1400 2    50   ~ 0
4
Text Label 12400 1500 2    50   ~ 0
5
Text Notes 12300 1200 2    50   ~ 0
Interline PWM 1
Text Label 6550 7350 2    50   ~ 0
2
Text Label 6550 8200 2    50   ~ 0
4
Text Label 6600 9050 2    50   ~ 0
6
Text Label 6600 9850 2    50   ~ 0
8
Text Label 7300 7350 0    50   ~ 0
3
Text Label 7300 8200 0    50   ~ 0
5
Text Label 7350 9050 0    50   ~ 0
7
Text Label 7350 9850 0    50   ~ 0
9
Text Notes 6450 7350 2    50   ~ 0
Interline PWM 1
Text Notes 6450 8200 2    50   ~ 0
Interline PWM 2
Text Notes 6500 9050 2    50   ~ 0
Interline PWM 3
Text Notes 6500 9850 2    50   ~ 0
Interline PWM 4
Text Notes 7400 7350 0    50   ~ 0
Analog select 1\n
Text Notes 7400 8200 0    50   ~ 0
Analog select 2\n
Text Notes 7450 9050 0    50   ~ 0
Analog select 3\n
Text Notes 7450 9850 0    50   ~ 0
Analog select 4\n
Wire Wire Line
	2650 10350 2650 10450
Wire Wire Line
	2650 10550 2650 10650
Wire Wire Line
	2650 10750 2650 10850
Wire Wire Line
	2650 10950 2650 11050
Text Notes 15000 1700 0    50   ~ 0
A/D IO 2
Text Notes 15000 2400 0    50   ~ 0
A/D IO 3
Text Notes 15350 3000 2    50   ~ 0
A/D IO 4
Text Notes 3400 10350 0    50   ~ 0
A/D IO 2
Text Notes 3450 11050 0    50   ~ 0
A/D IO 4
Text Label 2650 8700 0    50   ~ 0
21
Text Notes 2800 8700 0    50   ~ 0
A/D IO 1
Text Label 2650 8800 0    50   ~ 0
20
Text Notes 3150 8800 2    50   ~ 0
A/D IO 2
Text Label 2650 8900 0    50   ~ 0
19
Text Notes 2800 8900 0    50   ~ 0
A/D IO 3
Text Label 2650 9000 0    50   ~ 0
18
Text Notes 3150 9000 2    50   ~ 0
A/D IO 4
Text Notes 3900 8650 0    59   ~ 0
Distribute external analog inputs
Text Notes 15000 1300 0    50   ~ 0
LED pot 2
Text Notes 15000 1400 0    50   ~ 0
LED pot 3
Text Notes 15000 1500 0    50   ~ 0
LED pot 4
Text Notes 15000 1200 0    50   ~ 0
LED pot 1
Text Notes 11300 2600 0    50   ~ 0
Manual/Auto Switch
Text Notes 11400 2900 0    50   ~ 0
over-temp alarm
Text Notes 15000 1800 0    50   ~ 0
MOSFET temp 1
Text Notes 15000 1900 0    50   ~ 0
MOSFET temp 2
Text Notes 15000 2000 0    50   ~ 0
MOSFET temp 3
Text Notes 12050 3000 2    50   ~ 0
Ext fan PWM
Text Notes 12050 2000 2    50   ~ 0
Internal fan PWM
Text Notes 15000 2100 0    50   ~ 0
MOSFET temp 4
Text Notes 15000 2600 0    50   ~ 0
Resistor temp 1
Text Notes 15000 2700 0    50   ~ 0
Resistor temp 2
Text Notes 15000 2800 0    50   ~ 0
Resistor temp 3
Text Notes 15000 2900 0    50   ~ 0
Resistor temp 4
NoConn ~ 14700 2200
Text Notes 11400 2200 0    50   ~ 0
Status LED 4
Text Notes 11400 2100 0    50   ~ 0
Status LED 3
Text Label 12400 1900 2    50   ~ 0
9
Text Label 12400 1800 2    50   ~ 0
8
Text Label 12400 1700 2    50   ~ 0
7
Text Label 12400 1600 2    50   ~ 0
6
Text Notes 11350 1100 0    50   ~ 0
Status LED 2
Text Notes 11350 1000 0    50   ~ 0
Status LED 1
Text Notes 12300 1800 2    50   ~ 0
Interline PWM 4
Text Notes 12300 1600 2    50   ~ 0
Interline PWM 3
Text Notes 12300 1400 2    50   ~ 0
Interline PWM 2
Text Notes 11700 1300 0    50   ~ 0
Analog select 1\n
Text Notes 11700 1500 0    50   ~ 0
Analog select 2\n
Text Notes 11700 1700 0    50   ~ 0
Analog select 3\n
Text Notes 11700 1900 0    50   ~ 0
Analog select 4\n
$Comp
L Device:R_Pack08 RN?
U 1 1 5F57894B
P 2450 10750
F 0 "RN?" V 1833 10750 50  0000 C CNN
F 1 "4.7k" V 1924 10750 50  0000 C CNN
F 2 "Custom Footprints:Bourns_16-SOIC_8x_r-pack" V 2925 10750 50  0001 C CNN
F 3 "file://media.digikey.com/api/videos/videoplayer/smallplayer/N0706.pdf" H 2450 10750 50  0001 C CNN
F 4 "Bourns Inc." V 2450 10750 50  0001 C CNN "Manufacturer"
F 5 "4816P-T01-472LF" V 2450 10750 50  0001 C CNN "Part #"
	1    2450 10750
	0    1    1    0   
$EndComp
Text Notes 850  10150 1    59   ~ 0
I/O maximum voltage: 27V (160mW)\nNiDaq PCI-6110 is +/- 10V 5mA\nTherefore, minimum impedance is 2000 Ohms
Text Notes 4850 10650 0    50   ~ 0
A/D IO 3
Text Notes 4800 9950 0    50   ~ 0
A/D IO 1
Text Notes 3250 9700 0    59   ~ 0
4-channel analog/digital IO with 0-3.3V clamp
$Comp
L power:GND #PWR?
U 1 1 5F48EF1A
P 4700 11050
F 0 "#PWR?" H 4700 10800 50  0001 C CNN
F 1 "GND" V 4800 11000 50  0000 R CNN
F 2 "" H 4700 11050 50  0001 C CNN
F 3 "" H 4700 11050 50  0001 C CNN
	1    4700 11050
	0    -1   -1   0   
$EndComp
Text Label 4700 10850 0    50   ~ 0
3.3V
Text Label 3950 10850 2    50   ~ 0
3.3V
$Comp
L power:GND #PWR?
U 1 1 5F48EF12
P 3950 10650
F 0 "#PWR?" H 3950 10400 50  0001 C CNN
F 1 "GND" V 3850 10600 50  0000 R CNN
F 2 "" H 3950 10650 50  0001 C CNN
F 3 "" H 3950 10650 50  0001 C CNN
	1    3950 10650
	0    1    1    0   
$EndComp
$Comp
L Custom_parts:BAT54SDW D?
U 1 1 5F48EF0C
P 4150 10750
F 0 "D?" H 4325 11097 60  0000 C CNN
F 1 "BAT54SDW" H 4325 10991 60  0000 C CNN
F 2 "Package_TO_SOT_SMD:SOT-363_SC-70-6" H 4350 10950 60  0001 L CNN
F 3 "https://www.diodes.com/assets/Datasheets/ds11005.pdf" H 4350 11050 60  0001 L CNN
F 4 "Diodes Incorporated" H 4150 10750 50  0001 C CNN "Manufacturer"
F 5 "BAT54SDW-7-F" H 4150 10750 50  0001 C CNN "Part #"
	1    4150 10750
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F47759A
P 4650 10350
F 0 "#PWR?" H 4650 10100 50  0001 C CNN
F 1 "GND" V 4750 10300 50  0000 R CNN
F 2 "" H 4650 10350 50  0001 C CNN
F 3 "" H 4650 10350 50  0001 C CNN
	1    4650 10350
	0    -1   -1   0   
$EndComp
Text Label 4650 10150 0    50   ~ 0
3.3V
Text Label 3900 10150 2    50   ~ 0
3.3V
$Comp
L power:GND #PWR?
U 1 1 5F475E13
P 3900 9950
F 0 "#PWR?" H 3900 9700 50  0001 C CNN
F 1 "GND" V 3800 9900 50  0000 R CNN
F 2 "" H 3900 9950 50  0001 C CNN
F 3 "" H 3900 9950 50  0001 C CNN
	1    3900 9950
	0    1    1    0   
$EndComp
$Comp
L Custom_parts:BAT54SDW D?
U 1 1 5F46875C
P 4100 10050
F 0 "D?" H 4350 10350 60  0000 C CNN
F 1 "BAT54SDW" H 4350 10250 60  0000 C CNN
F 2 "Package_TO_SOT_SMD:SOT-363_SC-70-6" H 4300 10250 60  0001 L CNN
F 3 "https://www.diodes.com/assets/Datasheets/ds11005.pdf" H 4300 10350 60  0001 L CNN
F 4 "Diodes Incorporated" H 4100 10050 50  0001 C CNN "Manufacturer"
F 5 "BAT54SDW-7-F" H 4100 10050 50  0001 C CNN "Part #"
	1    4100 10050
	1    0    0    -1  
$EndComp
Text Notes 15000 1600 0    50   ~ 0
A/D IO 1
Text Notes 15850 3000 1    50   ~ 0
--ADC1--
Text Notes 15900 2400 2    50   ~ 0
ADC1
Text Notes 15900 2500 2    50   ~ 0
ADC0
Text Notes 15850 2100 1    50   ~ 0
-------ADC0------
Text Notes 15900 3100 2    50   ~ 0
ADC0
Text Notes 15900 3200 2    50   ~ 0
ADC0
Text Notes 11750 3100 2    50   ~ 0
ADC1
Text Notes 11750 3200 2    50   ~ 0
ADC1
Text Notes 15000 2500 0    50   ~ 0
Ext temp
Text Notes 16000 3600 1    50   ~ 0
https://forum.pjrc.com/attachment.php?attachmentid=10666&d=1495536536
$Comp
L Device:CP1_Small C?
U 1 1 5F41CA67
P 1150 3600
F 0 "C?" H 1100 3900 50  0000 L CNN
F 1 "14000uF" H 950 3800 50  0000 L CNN
F 2 "Capacitor_THT:CP_Radial_D18.0mm_P7.50mm" H 1150 3600 50  0001 C CNN
F 3 "https://media.digikey.com/pdf/Data%20Sheets/United%20Chemi-Con%20PDFs/LBK_Series.pdf" H 1150 3600 50  0001 C CNN
F 4 "United Chemi-Con" H 1150 3600 50  0001 C CNN "Manufacturer"
F 5 "ELBK250ELL143AM40S" H 1150 3600 50  0001 C CNN "Part #"
	1    1150 3600
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F41E64F
P 1150 3700
F 0 "#PWR?" H 1150 3450 50  0001 C CNN
F 1 "GND" H 1150 3700 50  0000 R CNN
F 2 "" H 1150 3700 50  0001 C CNN
F 3 "" H 1150 3700 50  0001 C CNN
	1    1150 3700
	1    0    0    -1  
$EndComp
$Comp
L Device:CP1_Small C?
U 1 1 5F422553
P 1150 4150
F 0 "C?" H 1100 4450 50  0000 L CNN
F 1 "14000uF" H 950 4350 50  0000 L CNN
F 2 "Capacitor_THT:CP_Radial_D18.0mm_P7.50mm" H 1150 4150 50  0001 C CNN
F 3 "https://media.digikey.com/pdf/Data%20Sheets/United%20Chemi-Con%20PDFs/LBK_Series.pdf" H 1150 4150 50  0001 C CNN
F 4 "United Chemi-Con" H 1150 4150 50  0001 C CNN "Manufacturer"
F 5 "ELBK250ELL143AM40S" H 1150 4150 50  0001 C CNN "Part #"
	1    1150 4150
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F422559
P 1150 4250
F 0 "#PWR?" H 1150 4000 50  0001 C CNN
F 1 "GND" H 1150 4250 50  0000 R CNN
F 2 "" H 1150 4250 50  0001 C CNN
F 3 "" H 1150 4250 50  0001 C CNN
	1    1150 4250
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F437161
P 1400 3600
F 0 "C?" H 1400 3300 50  0000 C CNN
F 1 "47uF" H 1400 3400 50  0000 C CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 1400 3600 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 1400 3600 50  0001 C CNN
F 4 "Murata Electronics" H 1400 3600 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 1400 3600 50  0001 C CNN "Part #"
	1    1400 3600
	-1   0    0    1   
$EndComp
$Comp
L Device:CP1_Small C?
U 1 1 5F45FB0D
P 1150 4700
F 0 "C?" H 1100 5000 50  0000 L CNN
F 1 "14000uF" H 950 4900 50  0000 L CNN
F 2 "Capacitor_THT:CP_Radial_D18.0mm_P7.50mm" H 1150 4700 50  0001 C CNN
F 3 "https://media.digikey.com/pdf/Data%20Sheets/United%20Chemi-Con%20PDFs/LBK_Series.pdf" H 1150 4700 50  0001 C CNN
F 4 "United Chemi-Con" H 1150 4700 50  0001 C CNN "Manufacturer"
F 5 "ELBK250ELL143AM40S" H 1150 4700 50  0001 C CNN "Part #"
	1    1150 4700
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F45FB13
P 1150 4800
F 0 "#PWR?" H 1150 4550 50  0001 C CNN
F 1 "GND" H 1150 4800 50  0000 R CNN
F 2 "" H 1150 4800 50  0001 C CNN
F 3 "" H 1150 4800 50  0001 C CNN
	1    1150 4800
	1    0    0    -1  
$EndComp
$Comp
L Device:CP1_Small C?
U 1 1 5F465A2E
P 1150 5250
F 0 "C?" H 1100 5550 50  0000 L CNN
F 1 "14000uF" H 950 5450 50  0000 L CNN
F 2 "Capacitor_THT:CP_Radial_D18.0mm_P7.50mm" H 1150 5250 50  0001 C CNN
F 3 "https://media.digikey.com/pdf/Data%20Sheets/United%20Chemi-Con%20PDFs/LBK_Series.pdf" H 1150 5250 50  0001 C CNN
F 4 "United Chemi-Con" H 1150 5250 50  0001 C CNN "Manufacturer"
F 5 "ELBK250ELL143AM40S" H 1150 5250 50  0001 C CNN "Part #"
	1    1150 5250
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F465A34
P 1150 5350
F 0 "#PWR?" H 1150 5100 50  0001 C CNN
F 1 "GND" H 1150 5350 50  0000 R CNN
F 2 "" H 1150 5350 50  0001 C CNN
F 3 "" H 1150 5350 50  0001 C CNN
	1    1150 5350
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F47CF6E
P 1400 3700
F 0 "#PWR?" H 1400 3450 50  0001 C CNN
F 1 "GND" H 1400 3700 50  0000 R CNN
F 2 "" H 1400 3700 50  0001 C CNN
F 3 "" H 1400 3700 50  0001 C CNN
	1    1400 3700
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F47E382
P 1400 4150
F 0 "C?" H 1400 3850 50  0000 C CNN
F 1 "47uF" H 1400 3950 50  0000 C CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 1400 4150 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 1400 4150 50  0001 C CNN
F 4 "Murata Electronics" H 1400 4150 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 1400 4150 50  0001 C CNN "Part #"
	1    1400 4150
	-1   0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F47E388
P 1400 4250
F 0 "#PWR?" H 1400 4000 50  0001 C CNN
F 1 "GND" H 1400 4250 50  0000 R CNN
F 2 "" H 1400 4250 50  0001 C CNN
F 3 "" H 1400 4250 50  0001 C CNN
	1    1400 4250
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F488538
P 1400 4700
F 0 "C?" H 1400 4400 50  0000 C CNN
F 1 "47uF" H 1400 4500 50  0000 C CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 1400 4700 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 1400 4700 50  0001 C CNN
F 4 "Murata Electronics" H 1400 4700 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 1400 4700 50  0001 C CNN "Part #"
	1    1400 4700
	-1   0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F48853E
P 1400 4800
F 0 "#PWR?" H 1400 4550 50  0001 C CNN
F 1 "GND" H 1400 4800 50  0000 R CNN
F 2 "" H 1400 4800 50  0001 C CNN
F 3 "" H 1400 4800 50  0001 C CNN
	1    1400 4800
	1    0    0    -1  
$EndComp
$Comp
L Device:C_Small C?
U 1 1 5F48D56B
P 1400 5250
F 0 "C?" H 1400 4950 50  0000 C CNN
F 1 "47uF" H 1400 5050 50  0000 C CNN
F 2 "Capacitor_SMD:C_1206_3216Metric" H 1400 5250 50  0001 C CNN
F 3 "https://www.digikey.com/en/ptm/m/murata-electronics-north-america/high-cap-multilayer-ceramic-capacitors?pn_sku=490-16268-1-ND&part_id=7363259" H 1400 5250 50  0001 C CNN
F 4 "Murata Electronics" H 1400 5250 50  0001 C CNN "Manufacturer"
F 5 "GRM31CR61E476ME44L" H 1400 5250 50  0001 C CNN "Part #"
	1    1400 5250
	-1   0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F48D571
P 1400 5350
F 0 "#PWR?" H 1400 5100 50  0001 C CNN
F 1 "GND" H 1400 5350 50  0000 R CNN
F 2 "" H 1400 5350 50  0001 C CNN
F 3 "" H 1400 5350 50  0001 C CNN
	1    1400 5350
	1    0    0    -1  
$EndComp
Wire Wire Line
	1400 3500 1150 3500
Wire Wire Line
	1400 4050 1150 4050
Wire Wire Line
	1400 4600 1150 4600
Wire Wire Line
	1400 5150 1150 5150
Wire Wire Line
	1150 5150 850  5150
Wire Wire Line
	850  5150 850  4600
Wire Wire Line
	850  3500 1150 3500
Connection ~ 1150 5150
Connection ~ 1150 3500
Wire Wire Line
	1150 4050 850  4050
Connection ~ 1150 4050
Connection ~ 850  4050
Wire Wire Line
	850  4050 850  3500
Wire Wire Line
	1150 4600 850  4600
Connection ~ 1150 4600
Connection ~ 850  4600
Wire Wire Line
	850  4600 850  4050
Wire Wire Line
	1200 1650 1400 1650
Wire Wire Line
	1400 2000 1700 2000
Wire Wire Line
	1400 2000 1400 1650
Wire Wire Line
	1400 1150 1050 1150
Connection ~ 1400 1650
Wire Wire Line
	1400 1650 1850 1650
Wire Wire Line
	1400 1650 1400 1150
$Comp
L Connector:RJ45_Shielded J?
U 1 1 5F4E6538
P 2650 3750
F 0 "J?" H 2707 4417 50  0000 C CNN
F 1 "RJ45_Shielded" H 2707 4326 50  0000 C CNN
F 2 "Custom Footprints:Molex_RJ45_855437001" V 2650 3775 50  0001 C CNN
F 3 "https://www.molex.com/pdm_docs/sd/855437001_sd.pdf" V 2650 3775 50  0001 C CNN
F 4 "Molex" H 2650 3750 50  0001 C CNN "Manufacturer"
F 5 "0855437001" H 2650 3750 50  0001 C CNN "Part #"
	1    2650 3750
	1    0    0    -1  
$EndComp
$Comp
L Device:Jumper_NO_Small JP?
U 1 1 5F4EBD0C
P 1850 3500
F 0 "JP?" H 1850 3600 50  0000 C CNN
F 1 "LED wire" H 1850 3400 50  0000 C CNN
F 2 "Connector_Wire:SolderWire-2.5sqmm_1x02_P7.2mm_D2.4mm_OD3.6mm" H 1850 3500 50  0001 C CNN
F 3 "~" H 1850 3500 50  0001 C CNN
	1    1850 3500
	1    0    0    -1  
$EndComp
Text Label 1500 3500 0    50   ~ 0
LED1+
Wire Wire Line
	1750 3500 1400 3500
Connection ~ 1400 3500
$Comp
L Device:Jumper_NO_Small JP?
U 1 1 5F4F518E
P 1850 4050
F 0 "JP?" H 1850 4150 50  0000 C CNN
F 1 "LED wire" H 1850 3950 50  0000 C CNN
F 2 "Connector_Wire:SolderWire-2.5sqmm_1x02_P7.2mm_D2.4mm_OD3.6mm" H 1850 4050 50  0001 C CNN
F 3 "~" H 1850 4050 50  0001 C CNN
	1    1850 4050
	1    0    0    -1  
$EndComp
Text Label 1500 4050 0    50   ~ 0
LED2+
Wire Wire Line
	1750 4050 1400 4050
$Comp
L Device:Jumper_NO_Small JP?
U 1 1 5F4FB4E9
P 1850 4600
F 0 "JP?" H 1850 4700 50  0000 C CNN
F 1 "LED wire" H 1850 4500 50  0000 C CNN
F 2 "Connector_Wire:SolderWire-2.5sqmm_1x02_P7.2mm_D2.4mm_OD3.6mm" H 1850 4600 50  0001 C CNN
F 3 "~" H 1850 4600 50  0001 C CNN
	1    1850 4600
	1    0    0    -1  
$EndComp
Text Label 1500 4600 0    50   ~ 0
LED3+
Wire Wire Line
	1750 4600 1400 4600
$Comp
L Device:Jumper_NO_Small JP?
U 1 1 5F501F26
P 1850 5150
F 0 "JP?" H 1850 5250 50  0000 C CNN
F 1 "LED wire" H 1850 5050 50  0000 C CNN
F 2 "Connector_Wire:SolderWire-2.5sqmm_1x02_P7.2mm_D2.4mm_OD3.6mm" H 1850 5150 50  0001 C CNN
F 3 "~" H 1850 5150 50  0001 C CNN
	1    1850 5150
	1    0    0    -1  
$EndComp
Text Label 1500 5150 0    50   ~ 0
LED4+
Wire Wire Line
	1750 5150 1400 5150
Text Label 1950 3500 0    50   ~ 0
LED1-
Text Label 1950 4050 0    50   ~ 0
LED2-
Text Label 1950 4600 0    50   ~ 0
LED3-
Text Label 1950 5150 0    50   ~ 0
LED4-
Text Label 3050 4050 0    50   ~ 0
LED1+
Text Label 3050 3950 0    50   ~ 0
LED1-
$Comp
L power:GND #PWR?
U 1 1 5F50FCB4
P 2650 4250
F 0 "#PWR?" H 2650 4000 50  0001 C CNN
F 1 "GND" H 2850 4250 50  0000 R CNN
F 2 "" H 2650 4250 50  0001 C CNN
F 3 "" H 2650 4250 50  0001 C CNN
	1    2650 4250
	1    0    0    -1  
$EndComp
Text Label 3050 5000 0    50   ~ 0
LED3+
Text Label 3050 5100 0    50   ~ 0
LED3-
$Comp
L power:GND #PWR?
U 1 1 5F523406
P 2650 5600
F 0 "#PWR?" H 2650 5350 50  0001 C CNN
F 1 "GND" H 2850 5600 50  0000 R CNN
F 2 "" H 2650 5600 50  0001 C CNN
F 3 "" H 2650 5600 50  0001 C CNN
	1    2650 5600
	1    0    0    -1  
$EndComp
Text Label 4150 3850 0    50   ~ 0
LED2+
Text Label 4150 3550 0    50   ~ 0
LED2-
$Comp
L power:GND #PWR?
U 1 1 5F52B343
P 3750 4250
F 0 "#PWR?" H 3750 4000 50  0001 C CNN
F 1 "GND" H 3950 4250 50  0000 R CNN
F 2 "" H 3750 4250 50  0001 C CNN
F 3 "" H 3750 4250 50  0001 C CNN
	1    3750 4250
	1    0    0    -1  
$EndComp
Text Label 4150 4800 0    50   ~ 0
LED4+
Text Label 4150 4700 0    50   ~ 0
LED4-
$Comp
L power:GND #PWR?
U 1 1 5F5331E6
P 3750 5600
F 0 "#PWR?" H 3750 5350 50  0001 C CNN
F 1 "GND" H 3950 5600 50  0000 R CNN
F 2 "" H 3750 5600 50  0001 C CNN
F 3 "" H 3750 5600 50  0001 C CNN
	1    3750 5600
	1    0    0    -1  
$EndComp
Text Label 3050 3850 0    50   ~ 0
LED1+
Text Label 3050 3650 0    50   ~ 0
LED1+
Text Label 3050 3450 0    50   ~ 0
LED1+
Text Label 3050 3750 0    50   ~ 0
LED1-
Text Label 3050 3550 0    50   ~ 0
LED1-
Text Label 3050 3350 0    50   ~ 0
LED1-
Text Label 4150 4050 0    50   ~ 0
LED2+
Text Label 4150 3650 0    50   ~ 0
LED2+
Text Label 4150 3450 0    50   ~ 0
LED2+
Text Label 4150 3950 0    50   ~ 0
LED2-
Text Label 4150 3750 0    50   ~ 0
LED2-
Text Label 4150 3350 0    50   ~ 0
LED2-
Text Label 3050 5400 0    50   ~ 0
LED3+
Text Label 3050 5200 0    50   ~ 0
LED3+
Text Label 3050 4800 0    50   ~ 0
LED3+
Text Label 3050 5300 0    50   ~ 0
LED3-
Text Label 3050 4900 0    50   ~ 0
LED3-
Text Label 3050 4700 0    50   ~ 0
LED3-
Text Label 4150 5400 0    50   ~ 0
LED4+
Text Label 4150 5200 0    50   ~ 0
LED4+
Text Label 4150 5000 0    50   ~ 0
LED4+
Text Label 4150 5300 0    50   ~ 0
LED4-
Text Label 4150 5100 0    50   ~ 0
LED4-
Text Label 4150 4900 0    50   ~ 0
LED4-
Connection ~ 1400 4050
Connection ~ 1400 4600
Connection ~ 1400 5150
Text Label 700  3500 2    50   ~ 0
Vin
Wire Wire Line
	700  3500 850  3500
Connection ~ 850  3500
$Comp
L Connector:RJ45_Shielded J?
U 1 1 5F5B6B8F
P 3750 3750
F 0 "J?" H 3807 4417 50  0000 C CNN
F 1 "RJ45_Shielded" H 3807 4326 50  0000 C CNN
F 2 "Custom Footprints:Molex_RJ45_855437001" V 3750 3775 50  0001 C CNN
F 3 "https://www.molex.com/pdm_docs/sd/855437001_sd.pdf" V 3750 3775 50  0001 C CNN
F 4 "Molex" H 3750 3750 50  0001 C CNN "Manufacturer"
F 5 "0855437001" H 3750 3750 50  0001 C CNN "Part #"
	1    3750 3750
	1    0    0    -1  
$EndComp
$Comp
L Connector:RJ45_Shielded J?
U 1 1 5F5B8657
P 2650 5100
F 0 "J?" H 2707 5767 50  0000 C CNN
F 1 "RJ45_Shielded" H 2707 5676 50  0000 C CNN
F 2 "Custom Footprints:Molex_RJ45_855437001" V 2650 5125 50  0001 C CNN
F 3 "https://www.molex.com/pdm_docs/sd/855437001_sd.pdf" V 2650 5125 50  0001 C CNN
F 4 "Molex" H 2650 5100 50  0001 C CNN "Manufacturer"
F 5 "0855437001" H 2650 5100 50  0001 C CNN "Part #"
	1    2650 5100
	1    0    0    -1  
$EndComp
$Comp
L Connector:RJ45_Shielded J?
U 1 1 5F5BA0DA
P 3750 5100
F 0 "J?" H 3807 5767 50  0000 C CNN
F 1 "RJ45_Shielded" H 3807 5676 50  0000 C CNN
F 2 "Custom Footprints:Molex_RJ45_855437001" V 3750 5125 50  0001 C CNN
F 3 "https://www.molex.com/pdm_docs/sd/855437001_sd.pdf" V 3750 5125 50  0001 C CNN
F 4 "Molex" H 3750 5100 50  0001 C CNN "Manufacturer"
F 5 "0855437001" H 3750 5100 50  0001 C CNN "Part #"
	1    3750 5100
	1    0    0    -1  
$EndComp
$Comp
L Connector:RJ45_Shielded J?
U 1 1 5F5E1891
P 1200 8050
F 0 "J?" H 1250 7300 50  0000 R CNN
F 1 "RJ45_Shielded" H 1450 7400 50  0000 R CNN
F 2 "Custom Footprints:Molex_RJ45_855437001" V 1200 8075 50  0001 C CNN
F 3 "https://www.molex.com/pdm_docs/sd/855437001_sd.pdf" V 1200 8075 50  0001 C CNN
F 4 "Molex" H 1200 8050 50  0001 C CNN "Manufacturer"
F 5 "0855437001" H 1200 8050 50  0001 C CNN "Part #"
	1    1200 8050
	1    0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F5F2844
P 1200 7550
F 0 "#PWR?" H 1200 7300 50  0001 C CNN
F 1 "GND" H 1450 7500 50  0000 R CNN
F 2 "" H 1200 7550 50  0001 C CNN
F 3 "" H 1200 7550 50  0001 C CNN
	1    1200 7550
	-1   0    0    1   
$EndComp
$Comp
L Connector:RJ45_Shielded J?
U 1 1 5F601BD3
P 1200 9350
F 0 "J?" H 1250 8600 50  0000 R CNN
F 1 "RJ45_Shielded" H 1450 8700 50  0000 R CNN
F 2 "Custom Footprints:Molex_RJ45_855437001" V 1200 9375 50  0001 C CNN
F 3 "https://www.molex.com/pdm_docs/sd/855437001_sd.pdf" V 1200 9375 50  0001 C CNN
F 4 "Molex" H 1200 9350 50  0001 C CNN "Manufacturer"
F 5 "0855437001" H 1200 9350 50  0001 C CNN "Part #"
	1    1200 9350
	1    0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F601BD9
P 1200 8850
F 0 "#PWR?" H 1200 8600 50  0001 C CNN
F 1 "GND" H 1450 8800 50  0000 R CNN
F 2 "" H 1200 8850 50  0001 C CNN
F 3 "" H 1200 8850 50  0001 C CNN
	1    1200 8850
	-1   0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F606C36
P 1600 7850
F 0 "#PWR?" H 1600 7600 50  0001 C CNN
F 1 "GND" V 1600 7700 50  0000 R CNN
F 2 "" H 1600 7850 50  0001 C CNN
F 3 "" H 1600 7850 50  0001 C CNN
	1    1600 7850
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F607DFB
P 1600 8050
F 0 "#PWR?" H 1600 7800 50  0001 C CNN
F 1 "GND" V 1600 7900 50  0000 R CNN
F 2 "" H 1600 8050 50  0001 C CNN
F 3 "" H 1600 8050 50  0001 C CNN
	1    1600 8050
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F60812B
P 1600 8250
F 0 "#PWR?" H 1600 8000 50  0001 C CNN
F 1 "GND" V 1600 8100 50  0000 R CNN
F 2 "" H 1600 8250 50  0001 C CNN
F 3 "" H 1600 8250 50  0001 C CNN
	1    1600 8250
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F6083B9
P 1600 8450
F 0 "#PWR?" H 1600 8200 50  0001 C CNN
F 1 "GND" V 1600 8300 50  0000 R CNN
F 2 "" H 1600 8450 50  0001 C CNN
F 3 "" H 1600 8450 50  0001 C CNN
	1    1600 8450
	0    -1   -1   0   
$EndComp
Wire Wire Line
	2250 8300 2250 7750
Wire Wire Line
	2250 7750 1600 7750
Wire Wire Line
	1600 7950 2200 7950
Wire Wire Line
	2200 7950 2200 8400
Wire Wire Line
	2200 8400 2250 8400
Wire Wire Line
	2250 8500 2150 8500
Wire Wire Line
	2150 8500 2150 8150
Wire Wire Line
	2150 8150 1600 8150
Wire Wire Line
	1600 8350 2100 8350
Wire Wire Line
	2100 8350 2100 8600
Wire Wire Line
	2100 8600 2250 8600
$Comp
L power:GND #PWR?
U 1 1 5F61E491
P 1600 9150
F 0 "#PWR?" H 1600 8900 50  0001 C CNN
F 1 "GND" V 1600 9000 50  0000 R CNN
F 2 "" H 1600 9150 50  0001 C CNN
F 3 "" H 1600 9150 50  0001 C CNN
	1    1600 9150
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F61EC78
P 1600 9350
F 0 "#PWR?" H 1600 9100 50  0001 C CNN
F 1 "GND" V 1600 9200 50  0000 R CNN
F 2 "" H 1600 9350 50  0001 C CNN
F 3 "" H 1600 9350 50  0001 C CNN
	1    1600 9350
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F61EEEF
P 1600 9550
F 0 "#PWR?" H 1600 9300 50  0001 C CNN
F 1 "GND" V 1600 9400 50  0000 R CNN
F 2 "" H 1600 9550 50  0001 C CNN
F 3 "" H 1600 9550 50  0001 C CNN
	1    1600 9550
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F61F25B
P 1600 9750
F 0 "#PWR?" H 1600 9500 50  0001 C CNN
F 1 "GND" V 1600 9600 50  0000 R CNN
F 2 "" H 1600 9750 50  0001 C CNN
F 3 "" H 1600 9750 50  0001 C CNN
	1    1600 9750
	0    -1   -1   0   
$EndComp
$Comp
L Connector:RJ45_Shielded J?
U 1 1 5F6566BD
P 1200 10650
F 0 "J?" H 1250 9900 50  0000 R CNN
F 1 "RJ45_Shielded" H 1450 10000 50  0000 R CNN
F 2 "Custom Footprints:Molex_RJ45_855437001" V 1200 10675 50  0001 C CNN
F 3 "https://www.molex.com/pdm_docs/sd/855437001_sd.pdf" V 1200 10675 50  0001 C CNN
F 4 "Molex" H 1200 10650 50  0001 C CNN "Manufacturer"
F 5 "0855437001" H 1200 10650 50  0001 C CNN "Part #"
	1    1200 10650
	1    0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F6566C3
P 1200 10150
F 0 "#PWR?" H 1200 9900 50  0001 C CNN
F 1 "GND" H 1450 10100 50  0000 R CNN
F 2 "" H 1200 10150 50  0001 C CNN
F 3 "" H 1200 10150 50  0001 C CNN
	1    1200 10150
	-1   0    0    1   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F6566C9
P 1600 10450
F 0 "#PWR?" H 1600 10200 50  0001 C CNN
F 1 "GND" V 1600 10300 50  0000 R CNN
F 2 "" H 1600 10450 50  0001 C CNN
F 3 "" H 1600 10450 50  0001 C CNN
	1    1600 10450
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F6566CF
P 1600 10650
F 0 "#PWR?" H 1600 10400 50  0001 C CNN
F 1 "GND" V 1600 10500 50  0000 R CNN
F 2 "" H 1600 10650 50  0001 C CNN
F 3 "" H 1600 10650 50  0001 C CNN
	1    1600 10650
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F6566D5
P 1600 10850
F 0 "#PWR?" H 1600 10600 50  0001 C CNN
F 1 "GND" V 1600 10700 50  0000 R CNN
F 2 "" H 1600 10850 50  0001 C CNN
F 3 "" H 1600 10850 50  0001 C CNN
	1    1600 10850
	0    -1   -1   0   
$EndComp
$Comp
L power:GND #PWR?
U 1 1 5F6566DB
P 1600 11050
F 0 "#PWR?" H 1600 10800 50  0001 C CNN
F 1 "GND" V 1600 10900 50  0000 R CNN
F 2 "" H 1600 11050 50  0001 C CNN
F 3 "" H 1600 11050 50  0001 C CNN
	1    1600 11050
	0    -1   -1   0   
$EndComp
Wire Wire Line
	1600 10350 2250 10350
Wire Wire Line
	1600 10550 2250 10550
Wire Wire Line
	1600 10750 2250 10750
Wire Wire Line
	1600 10950 2250 10950
Wire Wire Line
	2250 9650 2250 9000
Wire Wire Line
	1600 9650 2250 9650
Wire Wire Line
	2200 9450 2200 8900
Wire Wire Line
	2200 8900 2250 8900
Wire Wire Line
	1600 9450 2200 9450
Wire Wire Line
	2150 9250 2150 8800
Wire Wire Line
	2150 8800 2250 8800
Wire Wire Line
	1600 9250 2150 9250
Wire Wire Line
	2100 9050 2100 8700
Wire Wire Line
	2100 8700 2250 8700
Wire Wire Line
	1600 9050 2100 9050
Wire Wire Line
	3200 6800 3200 6550
Wire Wire Line
	3200 6550 3900 6550
Wire Wire Line
	2700 6800 3200 6800
Wire Wire Line
	3150 6900 3150 6250
Wire Wire Line
	2700 6900 3150 6900
Wire Wire Line
	3150 6250 5100 6250
$Comp
L Device:C_Small C?
U 1 1 5F6D2E0F
P 6450 7750
F 0 "C?" V 6650 7700 50  0000 L CNN
F 1 "2.2uF" V 6550 7600 50  0000 L CNN
F 2 "Custom Footprints:0508_Capacitor" H 6450 7750 50  0001 C CNN
F 3 "https://ds.yuden.co.jp/TYCOMPAS/ut/detail?pn=TWK212B7225MD-T%20&u=M" H 6450 7750 50  0001 C CNN
F 4 "Taiyo Yuden" H 6450 7750 50  0001 C CNN "Manufacturer"
F 5 "TWK212B7225MD-T" H 6450 7750 50  0001 C CNN "Part #"
	1    6450 7750
	0    1    1    0   
$EndComp
Connection ~ 6550 7750
$Comp
L power:GND #PWR?
U 1 1 5F6D4356
P 6350 7750
F 0 "#PWR?" H 6350 7500 50  0001 C CNN
F 1 "GND" V 6355 7622 50  0000 R CNN
F 2 "" H 6350 7750 50  0001 C CNN
F 3 "" H 6350 7750 50  0001 C CNN
	1    6350 7750
	0    1    1    0   
$EndComp
$EndSCHEMATC
