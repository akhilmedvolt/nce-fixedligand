from construction import fragment_construction
from scoring import fragment_scoring
from docking import fragment_docking

if __name__ == "__main__":
    fragment_scoring()
    print("Scoring Finished")
    print()
    print("------------------------------------------")
    fragment_construction()
    print("Construction Finished")
    print()
    print("------------------------------------------")
    fragment_docking()
    print("Execution Complete, Results Stored")
