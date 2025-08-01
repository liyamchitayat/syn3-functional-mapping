#!/usr/bin/env python3
"""
Download syn3A proteins - wrapper for data_downloader
"""

from data_downloader import DataDownloader

if __name__ == "__main__":
    downloader = DataDownloader("../../data")
    downloader.download_syn3a_proteins()