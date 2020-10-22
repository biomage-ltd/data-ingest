#!/usr/bin/python3
import os


def main():
    print("Biomage ingest running.")

    if os.listdir("/output"):
        print(os.listdir("/output"))
        print(
            "The output directory is not empty. Please clear it out before proceeding."
        )
        exit(1)

    print("Folders seem okay, starting R pre-processing...")


if __name__ == "__main__":
    main()