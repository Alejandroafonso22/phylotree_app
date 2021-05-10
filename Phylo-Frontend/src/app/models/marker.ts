/**
 * @file Model class to create the a marker object
 * @author Gerard Garcia
 * @version 1.0
 * @date 10/05/2021
*/
export class Marker {
    
    /**
     * Creates an instance of Marker.
     * @param {number} marker_id
     * @param {number} specie_id
     * @param {number} user_id
     * @param {number} longitude
     * @param {number} latitude
     * @param {Date} date
     * @param {string} time
     * @memberof Marker
     */
    constructor(private marker_id: number,
                private specie_id: number,
                private user_id: number,
                private longitude: number,
                private latitude: number,
                private date: Date,
                private time: string,
                ) {

    }
}