from genesys.assistants.ra import *
import genesys.DNAToolKit as bio_tools

import os

avaliable_functions = [
    {
        "type": "function",
        "function": {
            "name": "sequence_type",
            "description": "Get the type of sequences in a FASTA file.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    }
                },
                "required": ["filepath"]
            },
            "returns": {
                "type": "object",
                "properties": {
                    "sequence_type": {
                        "type": "string",
                        "enum": ["DNA", "RNA", "Protein"],
                    }
                },
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "count_occurrences",
            "description": "Count the number of nucleotides for each DNA/RNA sequence or amino acids for each protein in a FASTA file",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "transcription",
            "description": "Transcribe a DNA sequence to RNA.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "complementary",
            "description": "Find the complementary DNA sequence to a given DNA sequence.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "reverseComplementary",
            "description": "Find the reverse complementary of a sequence.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "gc_content",
            "description": "Calculate the GC content of a DNA/RNA sequence.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "translation",
            "description": "Translate a DNA sequence to a protein sequence.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "mass_calculator",
            "description": "Calculate the molecular mass of a DNA, RNA or protein sequence",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "restriction_sites",
            "description": "Find the restriction sites of a DNA sequence.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "isoelectric_point",
            "description": "Calculate the isoelectric point of a protein sequence.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "The protein sequence."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "multiple_sequence_alignment",
            "description": "Perform multiple sequence alignment using a FASTA file.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to the FASTA file."
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "construct_phylogenetic_tree",
            "description": "Construct a phylogenetic tree using a FASTA file.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "open_reading_frames",
            "description": "Finds the forward and reverse open reading frames",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "detect_snps",
            "description": "Detect single nucleotide polymorphisms (SNPs)",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                    },
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "find_motifs",
            "description": "Find motifs in a DNA sequence",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                    },
                    "motif": {
                        "type": "string",
                    }
                },
                "required": ["filepath", "motif"]
            }
        }
    }
]

assistant = ResearchAssistant(api_key=os.environ['OPENAI_API_KEY'])

user_input_medical = "Find me the restriction sites"
assistant.create_assistant_and_run(user_input=user_input_medical, assistant_name="Mendel",
                                   instructions="Call the respective functions to perform the bioinformatic operations",
                                   tools_list=avaliable_functions, thread_instructions="Your output should be the result of the function you called.",
                                   toolkit=bio_tools, file_path="/Users/rajagiddi/Code/DecoderAI/Bioinformatics/tests/fixtures/sequence.fasta")
